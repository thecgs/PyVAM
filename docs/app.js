const PYODIDE_URL = "https://cdn.jsdelivr.net/pyodide/v0.27.5/full/";
const PYVAM_FILES = ["config.py", "parserGB.py", "drawMT.py"];
const $ = (s) => document.querySelector(s);
const el = { files: $("#file-input"), selectedFiles: $("#selected-files"), chooseFiles: $("#choose-files"), accessions: $("#accession-input"), addAccessions: $("#add-accessions"), status: $("#runtime-status"), summary: $("#file-summary"), title: $("#plot-title"), note: $("#plot-note"), stats: $("#stats"), plot: $("#plot"), render: $("#render-button"), download: $("#download-button"), reset: $("#reset-button"), theme: $("#theme"), customColors: $("#custom-colors"), start: $("#start-feature"), labels: $("#show-labels") };
const state = { pyodide: null, files: [], inputs: [], nextInputId: 0, image: null, options: null, classOverrides: {}, pendingFiles: [], pendingAccessions: [] };
const mode = () => document.querySelector('input[name="view"]:checked').value;
const status = (text) => { el.status.textContent = text; };
const esc = (text) => String(text).replace(/[&<>"']/g, (c) => ({ "&":"&amp;", "<":"&lt;", ">":"&gt;", '"':"&quot;", "'":"&#039;" })[c]);
function hexTextColour(hex) { const value = hex.replace("#", ""); const brightness = (Number.parseInt(value.slice(0, 2), 16) * 299 + Number.parseInt(value.slice(2, 4), 16) * 587 + Number.parseInt(value.slice(4, 6), 16) * 114) / 1000; return brightness > 160 ? "#172a32" : "#ffffff"; }

const runtime = String.raw`
import base64, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pyvam import get_features, draw_circos_MT, draw_linear_MT, draw_linear_MT_nonproportional
from pyvam.config import MTColors, CommonNamesDict

GENE_CLASSES = {
    "Complex I (NADH dehydrogenase)": ["ND1", "ND2", "ND3", "ND4L", "ND4", "ND5", "ND6"],
    "Complex IV (Cytochrome c oxidase)": ["COX1", "COX2", "COX3"],
    "ATP synthase": ["ATPase6", "ATPase8"],
    "Cytochrome b": ["Cytb"],
    "transfer RNA": [name for name in next(iter(MTColors.values())) if name.startswith("tRNA-")],
    "ribosomal RNA": ["12S rRNA", "16S rRNA"],
    "D-loop (Non-coding region)": ["D-loop"],
    "Other genes": ["Other genes"],
}

def pyvam_options():
    return json.dumps({
        "themes": sorted(MTColors),
        "markers": sorted(set(CommonNamesDict.values()), key=str.casefold),
        "class_colours": {theme: {klass: palette.get(members[0], palette.get("Other genes", "#808080"))[:7] for klass, members in GENE_CLASSES.items()} for theme, palette in MTColors.items()},
    })

def pyvam_describe(paths_json):
    paths = json.loads(paths_json)
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

def pyvam_render(paths_json, view, theme, start, labels, custom_colors_json, options_json):
    # The image comes from PyVAM's public drawing API, not browser-side SVG code.
    paths = json.loads(paths_json)
    plt.close("all")
    colors = theme
    class_colours = json.loads(custom_colors_json)
    # Preserve each PyVAM theme exactly (including RGBA alpha) unless the user
    # explicitly changes a class colour in the web UI.
    colors = theme
    if class_colours:
        colors = dict(MTColors[theme.upper()])
        for klass, colour in class_colours.items():
            for gene in GENE_CLASSES[klass]: colors[gene] = colour
    web = json.loads(options_json)
    options = dict(colors=colors, start=start or None, isfilename2species=web["isfilename2species"],
                   abbr=web["abbr"], add_id=web["add_id"], remove_NCR=web["remove_ncr"], dpi=150)
    legend_position = web["legend_position"]
    if legend_position is not None: options["legend_postion"] = tuple(legend_position)
    if web["legend_size"] is not None: options["legend_size"] = web["legend_size"]
    if view == "circular":
        columns = min(3, len(paths)); rows = (len(paths) + columns - 1) // columns
        fig, axes = plt.subplots(rows, columns, subplot_kw={"projection": "polar"}, figsize=(6 * columns, 6 * rows), squeeze=False)
        for index, (path, axis) in enumerate(zip(paths, axes.flat)):
            circos_options = dict(options, show_gene_label=labels,
                                  show_legend=web["show_legend"] and index == len(paths) - 1,
                                  show_GC_circos=web["show_gc_circos"],
                                  gene_label_inner=web["gene_label_inner"],
                                  show_info=web["show_info"], direction=web["direction"],
                                  tidyname=web["tidyname"])
            if web["gene_label_size"] is not None: circos_options["gene_label_size"] = web["gene_label_size"]
            for key in ("radius", "info_fontsize", "GC_circos_height", "GC_circos_bin", "GC_circos_step"):
                if web[key] is not None: circos_options[key] = web[key]
            circos_options["GC_circos_color"] = web["GC_circos_color"]
            draw_circos_MT(path, axes=axis, **circos_options)
        for axis in list(axes.flat)[len(paths):]: axis.set_visible(False)
    elif view == "linear":
        linear_options = dict(options, show_gene_label=labels, show_legend=web["show_legend"], force_reoriented=web["force_reoriented"])
        for key in ("gene_label_size", "gene_label_color", "species_label_size", "species_label_color"):
            if web[key] is not None: linear_options[key] = web[key]
        fig, _ = draw_linear_MT(paths, **linear_options)
    else:
        order_options = dict(options, gene_label_size=web["gene_label_size"] if web["gene_label_size"] is not None else 9,
                             show_legend=web["show_legend"], force_reoriented=web["force_reoriented"])
        if not labels: order_options["gene_label_size"] = 0
        for key in ("gene_label_color", "species_label_size", "species_label_color", "height"):
            if web[key] is not None: order_options[key] = web[key]
        fig, _ = draw_linear_MT_nonproportional(paths, **order_options)
    output = "/home/pyvam-output.svg"
    fig.savefig(output, format="svg", bbox_inches="tight", dpi=150)
    plt.close(fig)
    with open(output, "rb") as handle: return base64.b64encode(handle.read()).decode("ascii")
`;

function colourEditors(theme) {
  state.classOverrides = {};
  const palette = state.options.class_colours[theme] || {};
  el.customColors.innerHTML = Object.entries(palette).map(([name, colour]) => `<label class="colour-row" style="--class-color:${esc(colour)}"><span class="class-name">${esc(name)}</span><input class="hex-colour" type="color" value="${esc(colour)}" data-class-colour="${esc(name)}" aria-label="${esc(name)} colour" /></label>`).join("");
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
    if (state.pendingFiles.length) await importFiles(state.pendingFiles.splice(0));
    if (state.pendingAccessions.length) await importAccessions(state.pendingAccessions.splice(0));
  } catch (error) {
    status("PyVAM runtime unavailable");
    el.summary.textContent = `Unable to load PyVAM: ${error.message}. This site must be published through the “Deploy PyVAM Web” GitHub Actions workflow, not directly from /docs.`;
    console.error(error);
  }
}

function updateStats(rows) {
  el.stats.innerHTML = rows.map((r) => `<div class="stat"><strong title="${esc(r.name)}">${esc(r.name)}</strong><span>${Number(r.length).toLocaleString()} bp · ${r.features} PyVAM features</span></div>`).join("");
}

function addInput(path, label) {
  state.inputs.push({ id: ++state.nextInputId, path, label });
  state.files = state.inputs.map((item) => item.path);
}

function renderSelectedFiles() {
  if (!state.inputs.length) { el.selectedFiles.textContent = "No GenBank files selected."; return; }
  el.selectedFiles.innerHTML = state.inputs.map((item) => `<span class="selected-file"><span title="${esc(item.label)}">${esc(item.label)}</span><button class="remove-file" type="button" data-input-id="${item.id}" aria-label="Remove ${esc(item.label)}">×</button></span>`).join("");
}

function optionalNumber(id) {
  const raw = $(id).value.trim();
  return raw === "" ? null : Number(raw);
}

function webOptions() {
  const x = optionalNumber("#legend-x"), y = optionalNumber("#legend-y");
  return {
    abbr: $("#abbr").checked,
    isfilename2species: $("#filename-species").checked,
    force_reoriented: $("#force-reoriented").checked,
    remove_ncr: $("#remove-ncr").checked,
    add_id: $("#add-id").checked,
    show_legend: $("#show-legend").checked,
    legend_size: optionalNumber("#legend-size"),
    legend_position: x === null && y === null ? null : [x === null ? 1 : x, y === null ? 0 : y],
    gene_label_size: optionalNumber("#gene-label-size"),
    gene_label_color: $("#gene-label-color").value,
    species_label_size: optionalNumber("#species-label-size"),
    species_label_color: $("#species-label-color").value,
    height: optionalNumber("#order-height"),
    radius: optionalNumber("#circos-radius"),
    gene_label_inner: $("#gene-label-inner").checked,
    show_info: $("#show-info").checked,
    info_fontsize: optionalNumber("#info-fontsize"),
    direction: Number($("#circos-direction").value),
    tidyname: $("#tidyname").checked,
    show_gc_circos: $("#show-gc-circos").checked,
    GC_circos_height: optionalNumber("#gc-circos-height"),
    GC_circos_color: $("#gc-circos-color").value,
    GC_circos_bin: optionalNumber("#gc-circos-bin"),
    GC_circos_step: optionalNumber("#gc-circos-step"),
  };
}

function updateViewOptions() {
  const view = mode();
  document.querySelectorAll(".linear-option").forEach((node) => { node.hidden = view === "circular"; });
  document.querySelectorAll(".order-option").forEach((node) => { node.hidden = view !== "order"; });
  document.querySelectorAll(".circular-option").forEach((node) => { node.hidden = view !== "circular"; });
  document.querySelectorAll(".gc-option").forEach((node) => { node.hidden = view !== "circular" || !$("#show-gc-circos").checked; });
}

async function importFiles(files) {
  if (!state.pyodide) { state.pendingFiles.push(...files); el.summary.textContent = `${state.pendingFiles.length} file(s) queued until PyVAM is ready.`; return; }
  for (const [i, file] of files.entries()) {
    // Preserve the original basename for isfilename2species; a unique parent
    // directory prevents collisions between separately selected same-name files.
    const uploadDir = `/home/uploads/${Date.now()}-${i}`;
    const name = file.name.replace(/[\\/]/g, "_") || "genome.gbk";
    state.pyodide.FS.mkdirTree(uploadDir);
    const path = `${uploadDir}/${name}`;
    state.pyodide.FS.writeFile(path, new Uint8Array(await file.arrayBuffer()));
    addInput(path, file.webkitRelativePath || file.name);
  }
  renderSelectedFiles();
  await refreshInputs();
}

async function importAccessions(accessions) {
  if (!state.pyodide) { state.pendingAccessions.push(...accessions); el.summary.textContent = `${state.pendingAccessions.length} accession(s) queued until PyVAM is ready.`; return; }
  for (const accession of accessions) {
    const clean = accession.trim(); if (!clean) continue;
    status(`Downloading ${clean} from NCBI…`);
    const endpoint = `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${encodeURIComponent(clean)}&rettype=gb&retmode=text`;
    const response = await fetch(endpoint);
    if (!response.ok) throw new Error(`NCBI returned HTTP ${response.status} for ${clean}`);
    state.pyodide.FS.writeFile(`/home/ncbi-${clean.replace(/[^A-Za-z0-9._-]+/g, "_")}.gb`, new Uint8Array(await response.arrayBuffer()));
    addInput(`/home/ncbi-${clean.replace(/[^A-Za-z0-9._-]+/g, "_")}.gb`, `NCBI: ${clean}`);
  }
  renderSelectedFiles();
  await refreshInputs();
}

async function refreshInputs() {
  try {
    state.pyodide.globals.set("web_paths_json", JSON.stringify(state.files));
    updateStats(JSON.parse(await state.pyodide.runPythonAsync("pyvam_describe(web_paths_json)")));
    el.summary.textContent = `${state.files.length} GenBank file${state.files.length === 1 ? "" : "s"} ready for PyVAM rendering.`;
    await render();
  } catch (error) { el.summary.textContent = `PyVAM could not parse the selected input(s): ${error.message}`; console.error(error); }
}

async function render() {
  if (!state.pyodide || !state.files.length) return;
  const hexInputs = [...document.querySelectorAll(".hex-colour")];
  if (!hexInputs.every(updateHexInput)) { el.summary.textContent = "Colours must use six-digit hexadecimal form, for example #FFEC00."; return; }
  const view = mode(), start = el.start.value;
  const customColors = state.classOverrides;
  el.render.disabled = true; status("PyVAM is rendering…");
  try {
    state.pyodide.globals.set("web_paths_json", JSON.stringify(state.files)); state.pyodide.globals.set("web_view", view); state.pyodide.globals.set("web_theme", el.theme.value); state.pyodide.globals.set("web_start", start); state.pyodide.globals.set("web_labels", el.labels.checked); state.pyodide.globals.set("web_custom_colors", JSON.stringify(customColors)); state.pyodide.globals.set("web_options", JSON.stringify(webOptions()));
    state.image = await state.pyodide.runPythonAsync("pyvam_render(web_paths_json, web_view, web_theme, web_start, web_labels, web_custom_colors, web_options)");
    el.plot.classList.remove("empty"); el.plot.innerHTML = `<img class="pyvam-figure" alt="PyVAM ${esc(view)} rendering" src="data:image/svg+xml;base64,${state.image}" />`;
    el.title.textContent = view === "order" ? "PyVAM gene-order comparison" : `PyVAM ${view} genome map${state.files.length === 1 ? "" : "s"}`;
    el.note.textContent = `${state.files.length} file${state.files.length === 1 ? "" : "s"} rendered by PyVAM + Matplotlib in this browser.`;
    el.download.disabled = false; status("PyVAM renderer ready");
  } catch (error) { el.summary.textContent = `PyVAM rendering failed: ${error.message}`; status("PyVAM render failed"); console.error(error); }
  finally { el.render.disabled = false; }
}

el.files.addEventListener("change", (event) => { const files = [...event.target.files]; event.target.value = ""; if (files.length) { el.selectedFiles.textContent = files.map((file) => file.webkitRelativePath || file.name).join("\n"); el.summary.textContent = `${files.length} file(s) selected; importing…`; } importFiles(files); });
el.selectedFiles.addEventListener("click", async (event) => {
  const button = event.target.closest("[data-input-id]"); if (!button) return;
  state.inputs = state.inputs.filter((item) => item.id !== Number(button.dataset.inputId));
  state.files = state.inputs.map((item) => item.path);
  renderSelectedFiles();
  if (state.files.length) await refreshInputs();
  else { state.image = null; el.summary.textContent = "No genomes loaded."; el.stats.innerHTML = ""; el.plot.className = "plot empty"; el.plot.innerHTML = '<div><span class="empty-icon">◎</span><p>Your PyVAM maps will appear here.</p></div>'; el.download.disabled = true; }
});
el.addAccessions.addEventListener("click", async () => {
  const accessions = el.accessions.value.split(/[\s,;]+/).filter(Boolean);
  if (!accessions.length) { el.summary.textContent = "Enter one or more NCBI nucleotide accession IDs first."; return; }
  try { await importAccessions(accessions); el.accessions.value = ""; }
  catch (error) { el.summary.textContent = `Could not retrieve NCBI accession(s): ${error.message}`; status("PyVAM renderer ready"); console.error(error); }
});
el.render.addEventListener("click", render);
document.querySelectorAll('input[name="view"]').forEach((input) => input.addEventListener("change", () => { updateViewOptions(); render(); }));
[el.start, el.labels, ...document.querySelectorAll(".option-group input, .option-group select")].forEach((input) => input.addEventListener("change", () => { updateViewOptions(); render(); }));
el.theme.addEventListener("change", () => { colourEditors(el.theme.value); render(); });
function updateHexInput(input) {
  const value = input.value.trim();
  if (!/^#[0-9a-fA-F]{6}$/.test(value)) { input.setCustomValidity("Use a six-digit hex colour, for example #FFEC00."); return false; }
  input.setCustomValidity(""); input.value = value.toUpperCase(); input.style.backgroundColor = input.value; input.style.color = hexTextColour(input.value);
  if (input.dataset.classColour) input.closest(".colour-row").style.setProperty("--class-color", input.value);
  return true;
}
el.customColors.addEventListener("input", (event) => { if (event.target.matches("[data-class-colour]") && updateHexInput(event.target)) { state.classOverrides[event.target.dataset.classColour] = event.target.value; render(); } });
el.customColors.addEventListener("change", (event) => { if (event.target.matches("[data-class-colour]") && updateHexInput(event.target)) { state.classOverrides[event.target.dataset.classColour] = event.target.value; render(); } });
document.querySelectorAll(".hex-colour").forEach((input) => { updateHexInput(input); input.addEventListener("input", () => { if (updateHexInput(input)) render(); }); });
el.reset.addEventListener("click", () => { state.files = []; state.inputs = []; state.pendingFiles = []; state.pendingAccessions = []; state.image = null; el.files.value = ""; renderSelectedFiles(); el.accessions.value = ""; el.summary.textContent = "No genomes loaded."; el.stats.innerHTML = ""; el.plot.className = "plot empty"; el.plot.innerHTML = '<div><span class="empty-icon">◎</span><p>Your PyVAM maps will appear here.</p></div>'; el.title.textContent = "Waiting for GenBank files"; el.note.textContent = "Load one or more annotated genomes to begin."; el.download.disabled = true; });
el.download.addEventListener("click", () => { if (!state.image) return; const link = Object.assign(document.createElement("a"), { href: `data:image/svg+xml;base64,${state.image}`, download: `pyvam-${mode()}-map.svg` }); link.click(); });
updateViewOptions();
initialise();
