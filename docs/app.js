const PYODIDE_URL = "https://cdn.jsdelivr.net/pyodide/v0.27.5/full/";
const PYVAM_FILES = ["config.py", "parserGB.py", "drawMT.py"];
const $ = (s) => document.querySelector(s);
const el = { files: $("#file-input"), status: $("#runtime-status"), summary: $("#file-summary"), title: $("#plot-title"), note: $("#plot-note"), stats: $("#stats"), plot: $("#plot"), render: $("#render-button"), download: $("#download-button"), reset: $("#reset-button"), theme: $("#theme"), customColors: $("#custom-colors"), start: $("#start-feature"), labels: $("#show-labels") };
const state = { pyodide: null, files: [], image: null, options: null };
const mode = () => document.querySelector('input[name="view"]:checked').value;
const status = (text) => { el.status.textContent = text; };
const esc = (text) => String(text).replace(/[&<>"']/g, (c) => ({ "&":"&amp;", "<":"&lt;", ">":"&gt;", '"':"&quot;", "'":"&#039;" })[c]);

const runtime = String.raw`
import base64, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pyvam import get_features, draw_circos_MT, draw_linear_MT, draw_linear_MT_nonproportional
from pyvam.config import MTColors, CommonNamesDict

def pyvam_options():
    return json.dumps({
        "themes": sorted(MTColors),
        "markers": sorted(set(CommonNamesDict.values()), key=str.casefold),
        "palettes": MTColors,
    })

def pyvam_describe(paths):
    result = []
    for path in paths:
        features = get_features(path, isfilename2species=True)
        source = features[0] if features else None
        result.append({
            "name": source.name if source is not None else os.path.basename(path),
            "length": len(source.location) if source is not None else 0,
            "features": max(len(features) - 1, 0),
        })
    return json.dumps(result)

def pyvam_render(paths, view, theme, start, labels, custom_colors_json):
    # The image comes from PyVAM's public drawing API, not browser-side SVG code.
    plt.close("all")
    colors = theme
    custom_colors = json.loads(custom_colors_json)
    if custom_colors:
        colors = dict(MTColors[theme.upper()])
        colors.update(custom_colors)
    options = dict(colors=colors, start=start or None, isfilename2species=True, dpi=150)
    if view == "circular":
        columns = min(3, len(paths)); rows = (len(paths) + columns - 1) // columns
        fig, axes = plt.subplots(rows, columns, subplot_kw={"projection": "polar"}, figsize=(6 * columns, 6 * rows), squeeze=False)
        for path, axis in zip(paths, axes.flat):
            draw_circos_MT(path, axes=axis, show_gene_label=labels, show_legend=False, show_GC_circos=True, **options)
        for axis in list(axes.flat)[len(paths):]: axis.set_visible(False)
    elif view == "linear":
        fig, _ = draw_linear_MT(paths, show_gene_label=labels, show_legend=True, force_reoriented=bool(start), **options)
    else:
        fig, _ = draw_linear_MT_nonproportional(paths, gene_label_size=9 if labels else 0, show_legend=True, force_reoriented=bool(start), **options)
    output = "/home/pyvam-output.svg"
    fig.savefig(output, format="svg", bbox_inches="tight", dpi=150)
    plt.close(fig)
    with open(output, "rb") as handle: return base64.b64encode(handle.read()).decode("ascii")
`;

function colourEditors(theme) {
  const palette = state.options.palettes[theme] || {};
  el.customColors.innerHTML = Object.entries(palette).filter(([name]) => name !== "source").map(([name, colour]) => `<label class="colour-row"><input type="checkbox" data-colour-enabled="${esc(name)}" /><span>${esc(name)}</span><input type="color" value="${esc(colour)}" data-colour-value="${esc(name)}" aria-label="${esc(name)} colour" /></label>`).join("");
}

function populatePyvamOptions(options) {
  state.options = options;
  el.theme.innerHTML = options.themes.map((name) => `<option value="${esc(name)}">${esc(name)}</option>`).join("");
  el.theme.value = options.themes.includes("MITOFISH") ? "MITOFISH" : options.themes[0];
  el.theme.disabled = false;
  el.start.innerHTML = `<option value="">GenBank origin</option>${options.markers.map((name) => `<option value="${esc(name)}">${esc(name)}</option>`).join("")}`;
  el.start.disabled = false;
  colourEditors(el.theme.value);
}

async function initialise() {
  try {
    await new Promise((ok, fail) => { const tag = document.createElement("script"); tag.src = `${PYODIDE_URL}pyodide.js`; tag.onload = ok; tag.onerror = fail; document.head.append(tag); });
    state.pyodide = await window.loadPyodide({ indexURL: PYODIDE_URL });
    status("Loading PyVAM runtime…");
    await state.pyodide.loadPackage(["matplotlib", "biopython"]);
    state.pyodide.FS.mkdirTree("/home/pyvam");
    for (const name of PYVAM_FILES) {
      const response = await fetch(`./pyvam/${name}`);
      if (!response.ok) throw new Error(`bundled PyVAM source unavailable: ${name}`);
      state.pyodide.FS.writeFile(`/home/pyvam/${name}`, await response.text());
    }
    state.pyodide.FS.writeFile("/home/pyvam/__init__.py", "from .parserGB import get_features\nfrom .drawMT import draw_circos_MT, draw_linear_MT, draw_linear_MT_nonproportional\n");
    await state.pyodide.runPythonAsync("import sys; sys.path.insert(0, '/home')");
    await state.pyodide.runPythonAsync(runtime);
    populatePyvamOptions(JSON.parse(await state.pyodide.runPythonAsync("pyvam_options()")));
    status("PyVAM renderer ready");
  } catch (error) {
    status("PyVAM runtime unavailable");
    el.summary.textContent = `Unable to load PyVAM: ${error.message}. This site must be published through the “Deploy PyVAM Web” GitHub Actions workflow, not directly from /docs.`;
    console.error(error);
  }
}

function updateStats(rows) {
  el.stats.innerHTML = rows.map((r) => `<div class="stat"><strong title="${esc(r.name)}">${esc(r.name)}</strong><span>${Number(r.length).toLocaleString()} bp · ${r.features} PyVAM features</span></div>`).join("");
}

async function importFiles(files) {
  if (!state.pyodide) { el.summary.textContent = "PyVAM is still loading; please try again in a moment."; return; }
  state.files = [];
  for (const [i, file] of files.entries()) {
    const name = `${i}-${file.name.replace(/[^A-Za-z0-9._-]+/g, "_") || "genome.gbk"}`;
    state.pyodide.FS.writeFile(`/home/${name}`, await file.arrayBuffer());
    state.files.push(`/home/${name}`);
  }
  try {
    state.pyodide.globals.set("web_paths", state.files);
    updateStats(JSON.parse(await state.pyodide.runPythonAsync("pyvam_describe(web_paths)")));
    el.summary.textContent = `${state.files.length} GenBank file${state.files.length === 1 ? "" : "s"} ready for PyVAM rendering.`;
    await render();
  } catch (error) { el.summary.textContent = `PyVAM could not parse the selected file(s): ${error.message}`; console.error(error); }
}

async function render() {
  if (!state.pyodide || !state.files.length) return;
  const view = mode(), start = el.start.value;
  const customColors = Object.fromEntries([...el.customColors.querySelectorAll("[data-colour-enabled]")]
    .filter((box) => box.checked)
    .map((box) => [box.dataset.colourEnabled, el.customColors.querySelector(`[data-colour-value="${CSS.escape(box.dataset.colourEnabled)}"]`).value]));
  el.render.disabled = true; status("PyVAM is rendering…");
  try {
    state.pyodide.globals.set("web_paths", state.files); state.pyodide.globals.set("web_view", view); state.pyodide.globals.set("web_theme", el.theme.value); state.pyodide.globals.set("web_start", start); state.pyodide.globals.set("web_labels", el.labels.checked); state.pyodide.globals.set("web_custom_colors", JSON.stringify(customColors));
    state.image = await state.pyodide.runPythonAsync("pyvam_render(web_paths, web_view, web_theme, web_start, web_labels, web_custom_colors)");
    el.plot.classList.remove("empty"); el.plot.innerHTML = `<img class="pyvam-figure" alt="PyVAM ${esc(view)} rendering" src="data:image/svg+xml;base64,${state.image}" />`;
    el.title.textContent = view === "order" ? "PyVAM gene-order comparison" : `PyVAM ${view} genome map${state.files.length === 1 ? "" : "s"}`;
    el.note.textContent = `${state.files.length} file${state.files.length === 1 ? "" : "s"} rendered by PyVAM + Matplotlib in this browser.`;
    el.download.disabled = false; status("PyVAM renderer ready");
  } catch (error) { el.summary.textContent = `PyVAM rendering failed: ${error.message}`; status("PyVAM render failed"); console.error(error); }
  finally { el.render.disabled = false; }
}

el.files.addEventListener("change", (event) => importFiles([...event.target.files]));
el.render.addEventListener("click", render);
[...document.querySelectorAll('input[name="view"]'), el.start, el.labels].forEach((input) => input.addEventListener("change", render));
el.theme.addEventListener("change", () => { colourEditors(el.theme.value); render(); });
el.customColors.addEventListener("input", render);
el.customColors.addEventListener("change", render);
el.reset.addEventListener("click", () => { state.files = []; state.image = null; el.files.value = ""; el.summary.textContent = "No genomes loaded."; el.stats.innerHTML = ""; el.plot.className = "plot empty"; el.plot.innerHTML = '<div><span class="empty-icon">◎</span><p>Your PyVAM maps will appear here.</p></div>'; el.title.textContent = "Waiting for GenBank files"; el.note.textContent = "Load one or more annotated genomes to begin."; el.download.disabled = true; });
el.download.addEventListener("click", () => { if (!state.image) return; const link = Object.assign(document.createElement("a"), { href: `data:image/svg+xml;base64,${state.image}`, download: `pyvam-${mode()}-map.svg` }); link.click(); });
initialise();
