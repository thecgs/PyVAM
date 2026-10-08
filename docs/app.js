const PYODIDE_URL = "https://cdn.jsdelivr.net/pyodide/v0.27.5/full/";

const parserCode = String.raw`
import json, re

CODING = {
    "COX1", "COX2", "COX3", "CYTB", "ATP6", "ATP8", "ND1", "ND2", "ND3",
    "ND4", "ND4L", "ND5", "ND6"
}

def _canon(raw, feature_type):
    value = (raw or "").strip().upper()
    compact = re.sub(r"[^A-Z0-9]", "", value)
    if "CONTROL REGION" in value or "D-LOOP" in value or "D LOOP" in value:
        return "D-loop"
    if "12S" in value or "RNR1" in value: return "12S rRNA"
    if "16S" in value or "RNR2" in value: return "16S rRNA"
    if "TRNA" in value:
        amino = {"ALANINE":"Ala", "CYSTEINE":"Cys", "ASPARTIC":"Asp", "GLUTAMIC":"Glu", "PHENYLALANINE":"Phe", "GLYCINE":"Gly", "HISTIDINE":"His", "ISOLEUCINE":"Ile", "LYSINE":"Lys", "LEUCINE":"Leu", "METHIONINE":"Met", "ASPARAGINE":"Asn", "PROLINE":"Pro", "GLUTAMINE":"Gln", "ARGININE":"Arg", "SERINE":"Ser", "THREONINE":"Thr", "VALINE":"Val", "TRYPTOPHAN":"Trp", "TYROSINE":"Tyr"}
        for full_name, short_name in amino.items():
            if full_name in compact: return "tRNA-" + short_name
        three_letter = {"ALA":"Ala", "CYS":"Cys", "ASP":"Asp", "GLU":"Glu", "PHE":"Phe", "GLY":"Gly", "HIS":"His", "ILE":"Ile", "LYS":"Lys", "LEU":"Leu", "MET":"Met", "ASN":"Asn", "PRO":"Pro", "GLN":"Gln", "ARG":"Arg", "SER":"Ser", "THR":"Thr", "VAL":"Val", "TRP":"Trp", "TYR":"Tyr"}
        if compact[4:] in three_letter: return "tRNA-" + three_letter[compact[4:]]
        code = re.search(r"TRNA[^A-Z]*([ACDEFGHIKLMNPQRSTVWY])$", value)
        one_letter = {"A":"Ala", "C":"Cys", "D":"Asp", "E":"Glu", "F":"Phe", "G":"Gly", "H":"His", "I":"Ile", "K":"Lys", "L":"Leu", "M":"Met", "N":"Asn", "P":"Pro", "Q":"Gln", "R":"Arg", "S":"Ser", "T":"Thr", "V":"Val", "W":"Trp", "Y":"Tyr"}
        return "tRNA-" + one_letter.get(code.group(1), "?") if code else "tRNA"
    aliases = [("COX", "COX"), ("COI", "COX1"), ("COII", "COX2"), ("COIII", "COX3"), ("CYTB", "Cytb"), ("ATPASE", "ATP"), ("ATP", "ATP"), ("NADH", "ND"), ("ND", "ND")]
    for prefix, target in aliases:
        if prefix in compact:
            match = re.search(r"(?:SUBUNIT)?(\dL?|L)", compact[compact.find(prefix):])
            if match:
                number = match.group(1)
                if target == "COX": return "COX" + number
                if target == "ATP": return "ATPase" + number
                return "ND" + number
    if feature_type.lower() in ("d_loop", "control_region"): return "D-loop"
    return (raw or feature_type).strip() or feature_type

def _location(location):
    pairs = [(int(a) - 1, int(b)) for a, b in re.findall(r"(\d+)\.\.(\d+)", location)]
    if not pairs: return None
    # Keep join() parts separately: a feature spanning the circular origin must
    # remain two segments rather than becoming one near-genome-length interval.
    return min(a for a, b in pairs), max(b for a, b in pairs), -1 if "complement" in location.lower() else 1, pairs

def parse_genbank(text):
    genomes = []
    for block in re.split(r"\n//\s*", text):
        if "FEATURES" not in block: continue
        header, feature_text = block.split("FEATURES", 1)
        locus = re.search(r"^LOCUS\s+(\S+).*?(\d+)\s+bp", header, re.M)
        sequence_match = re.search(r"^ORIGIN\s*\n([\s\S]*)$", feature_text, re.M)
        sequence = re.sub(r"[^A-Za-z]", "", sequence_match.group(1)).upper() if sequence_match else ""
        organism = re.search(r"^\s{2}ORGANISM\s+(.+)$", header, re.M)
        name = organism.group(1).strip() if organism else (locus.group(1) if locus else "Untitled genome")
        length = int(locus.group(2)) if locus else len(sequence)
        features, current = [], None
        for line in feature_text.splitlines():
            if line.startswith("ORIGIN"): break
            m = re.match(r"^\s{5}(\S+)\s+(.+)$", line)
            if m:
                if current: features.append(current)
                current = {"rawType":m.group(1), "location":m.group(2), "qualifiers":{}}
                continue
            q = re.match(r'^\s{21}/([^=]+)(?:=(.*))?$', line)
            if q and current:
                key, value = q.group(1), (q.group(2) or "").strip().strip('"')
                current["qualifiers"].setdefault(key, []).append(value)
            elif current and line.startswith("                     ") and current["qualifiers"]:
                key = list(current["qualifiers"])[-1]
                current["qualifiers"][key][-1] += line.strip().strip('"')
        if current: features.append(current)
        selected, seen = [], set()
        for feature in features:
            raw_type = feature["rawType"]
            if raw_type == "source": continue
            loc = _location(feature["location"])
            if not loc: continue
            values = feature["qualifiers"]
            raw_name = next((values[key][0] for key in ("gene", "product", "label", "note") if values.get(key)), raw_type)
            label = _canon(raw_name, raw_type)
            ftype = raw_type.lower()
            if ftype == "cds": kind = "cds"
            elif ftype == "trna" or label.startswith("tRNA-"): kind = "trna"
            elif ftype == "rrna" or "rRNA" in label: kind = "rrna"
            elif label == "D-loop" or ftype in ("d_loop", "control_region"): kind = "control"
            else: kind = "other"
            key = (label, loc[0], loc[1])
            # A GenBank annotation often has both gene and CDS entries for one locus.
            if key in seen: continue
            seen.add(key)
            selected.append({"label":label, "start":loc[0], "end":loc[1], "strand":loc[2], "parts":loc[3], "kind":kind, "rawType":raw_type})
        genomes.append({"name":name, "length":length, "sequenceLength":len(sequence), "features":selected})
    return json.dumps(genomes)
`;

const $ = (selector) => document.querySelector(selector);
const elements = {
  files: $("#file-input"), status: $("#runtime-status"), summary: $("#file-summary"),
  title: $("#plot-title"), note: $("#plot-note"), stats: $("#stats"), plot: $("#plot"),
  render: $("#render-button"), download: $("#download-button"), reset: $("#reset-button"),
  theme: $("#theme"), start: $("#start-feature"), labels: $("#show-labels"), strands: $("#show-strands"),
};
const state = { pyodide: null, genomes: [] };
const themes = {
  classic: { cds: "#e85d3f", trna: "#3e8d94", rrna: "#885f9c", control: "#e6b64c", other: "#9aa9a0" },
  muted: { cds: "#b56f5d", trna: "#668d83", rrna: "#89738f", control: "#c5a961", other: "#98a09b" },
  contrast: { cds: "#e63946", trna: "#0077b6", rrna: "#7209b7", control: "#f4a261", other: "#6c757d" },
};

function esc(text) { return String(text).replace(/[&<>"']/g, c => ({ "&":"&amp;", "<":"&lt;", ">":"&gt;", '"':"&quot;", "'":"&#039;" })[c]); }
function viewMode() { return document.querySelector('input[name="view"]:checked').value; }
function validFeatures(genome) { return genome.features.filter(f => f.end > f.start && f.start >= 0 && f.end <= genome.length); }
function rotatePosition(position, genome, start) {
  if (start === "origin") return position;
  const feature = genome.features.find(f => f.label.toUpperCase() === start.toUpperCase());
  return feature ? (position - feature.start + genome.length) % genome.length : position;
}
function segments(feature, genome, start) {
  return (feature.parts || [[feature.start, feature.end]]).flatMap(([from, to]) => {
    const a = rotatePosition(from, genome, start), b = rotatePosition(to, genome, start);
    return b > a ? [[a, b]] : [[a, genome.length], [0, b]];
  });
}
function polar(cx, cy, radius, angle) { return [cx + radius * Math.cos(angle), cy + radius * Math.sin(angle)]; }
function arcPath(cx, cy, inner, outer, start, end) {
  const large = end - start > Math.PI ? 1 : 0;
  const [a, b] = polar(cx, cy, outer, start), [c, d] = polar(cx, cy, outer, end), [e, f] = polar(cx, cy, inner, end), [g, h] = polar(cx, cy, inner, start);
  return `M ${a} ${b} A ${outer} ${outer} 0 ${large} 1 ${c} ${d} L ${e} ${f} A ${inner} ${inner} 0 ${large} 0 ${g} ${h} Z`;
}
function featureFill(feature, palette) { return palette[feature.kind] || palette.other; }

function circularSVG(genomes, palette, start, labels, showStrands) {
  const columns = Math.min(3, genomes.length), size = 340, rows = Math.ceil(genomes.length / columns);
  let svg = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${columns * size} ${rows * size}" role="img" aria-label="Circular mitochondrial genome maps"><style>.glabel{font:600 11px system-ui;fill:#1d3034}.gname{font:700 15px system-ui;fill:#172f35}.sub{font:12px system-ui;fill:#60727a}.reverse{stroke:#142f35;stroke-width:1.2}</style>`;
  genomes.forEach((genome, index) => {
    const cx = (index % columns) * size + size / 2, cy = Math.floor(index / columns) * size + 145, radius = 92;
    svg += `<g><circle cx="${cx}" cy="${cy}" r="${radius}" fill="none" stroke="#e1e7df" stroke-width="29"/>`;
    validFeatures(genome).forEach(feature => segments(feature, genome, start).forEach(([from, to]) => {
      const a = from / genome.length * 2 * Math.PI - Math.PI / 2, b = to / genome.length * 2 * Math.PI - Math.PI / 2;
      svg += `<path d="${arcPath(cx, cy, radius - 14, radius + 14, a, b)}" fill="${featureFill(feature, palette)}" ${showStrands && feature.strand < 0 ? 'class="reverse"' : ""}/>`;
    }));
    if (labels) validFeatures(genome).forEach(feature => {
      const mid = rotatePosition((feature.start + feature.end) / 2, genome, start) / genome.length * 2 * Math.PI - Math.PI / 2;
      const [x, y] = polar(cx, cy, radius + 34, mid);
      const anchor = Math.cos(mid) > .25 ? "start" : Math.cos(mid) < -.25 ? "end" : "middle";
      svg += `<text class="glabel" x="${x}" y="${y}" text-anchor="${anchor}" dominant-baseline="middle">${esc(feature.label)}</text>`;
    });
    svg += `<text class="gname" x="${cx}" y="${cy - 5}" text-anchor="middle">${esc(genome.name).slice(0, 28)}</text><text class="sub" x="${cx}" y="${cy + 15}" text-anchor="middle">${genome.length.toLocaleString()} bp · ${genome.features.length} features</text></g>`;
  });
  return svg + "</svg>";
}

function linearSVG(genomes, palette, start, labels, showStrands) {
  const width = 1060, row = labels ? 108 : 76, height = genomes.length * row + 55;
  let svg = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" role="img" aria-label="Linear mitochondrial genome maps"><style>.glabel{font:600 11px system-ui;fill:#21383c}.gname{font:700 14px system-ui;fill:#172f35}.axis{stroke:#cfdad2;stroke-width:2}.reverse{stroke:#142f35;stroke-width:1.2}</style>`;
  genomes.forEach((genome, i) => {
    const y = 38 + i * row, x = 180, length = 820;
    svg += `<text class="gname" x="12" y="${y + 5}">${esc(genome.name).slice(0, 24)}</text><text x="12" y="${y + 22}" fill="#60727a" font-size="11">${genome.length.toLocaleString()} bp</text><line class="axis" x1="${x}" x2="${x + length}" y1="${y}" y2="${y}"/>`;
    validFeatures(genome).forEach(feature => segments(feature, genome, start).forEach(([from, to]) => {
      const fx = x + from / genome.length * length, fw = Math.max(2, (to - from) / genome.length * length);
      svg += `<rect x="${fx}" y="${y - 11}" width="${fw}" height="22" rx="2" fill="${featureFill(feature, palette)}" ${showStrands && feature.strand < 0 ? 'class="reverse"' : ""}/>`;
    }));
    if (labels) validFeatures(genome).forEach(feature => { const pos = rotatePosition((feature.start + feature.end) / 2, genome, start); svg += `<text class="glabel" x="${x + pos / genome.length * length}" y="${y + 30}" text-anchor="middle">${esc(feature.label)}</text>`; });
  });
  return svg + "</svg>";
}

function orderSVG(genomes, palette, labels, showStrands) {
  const width = 1060, row = 83, height = genomes.length * row + 45;
  let svg = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" role="img" aria-label="Gene order comparison"><style>.glabel{font:600 10px system-ui;fill:#fff}.gname{font:700 14px system-ui;fill:#172f35}.reverse{stroke:#142f35;stroke-width:1.3}</style>`;
  genomes.forEach((genome, i) => {
    const y = 30 + i * row, features = validFeatures(genome).sort((a,b) => a.start-b.start), x = 190, available = 805;
    svg += `<text class="gname" x="12" y="${y + 6}">${esc(genome.name).slice(0, 24)}</text><text x="12" y="${y + 24}" fill="#60727a" font-size="11">${features.length} annotated loci</text>`;
    features.forEach((feature, n) => { const w = available / Math.max(features.length, 1), fx = x + n * w; svg += `<rect x="${fx}" y="${y - 17}" width="${Math.max(w - 1, 1)}" height="35" rx="2" fill="${featureFill(feature,palette)}" ${showStrands && feature.strand < 0 ? 'class="reverse"' : ""}/>${labels && w > 26 ? `<text class="glabel" x="${fx+w/2}" y="${y+5}" text-anchor="middle">${esc(feature.label).slice(0,10)}</text>` : ""}`; });
  });
  return svg + "</svg>";
}

function updateStats() {
  elements.stats.innerHTML = state.genomes.map(g => `<div class="stat"><strong title="${esc(g.name)}">${esc(g.name)}</strong><span>${g.length.toLocaleString()} bp · ${g.features.length} mapped features</span></div>`).join("");
}
function render() {
  if (!state.genomes.length) return;
  const mode = viewMode(), palette = themes[elements.theme.value], start = elements.start.value, labels = elements.labels.checked, strands = elements.strands.checked;
  elements.plot.classList.remove("empty");
  elements.plot.innerHTML = mode === "circular" ? circularSVG(state.genomes, palette, start, labels, strands) : mode === "linear" ? linearSVG(state.genomes, palette, start, labels, strands) : orderSVG(state.genomes, palette, labels, strands);
  elements.title.textContent = mode === "order" ? "Gene-order comparison" : `${mode[0].toUpperCase() + mode.slice(1)} genome maps`;
  elements.note.textContent = `${state.genomes.length} genome${state.genomes.length === 1 ? "" : "s"} · ${mode} presentation`;
  elements.download.disabled = false;
  updateStats();
}

async function loadPyodide() {
  try {
    await new Promise((resolve, reject) => { const script = document.createElement("script"); script.src = `${PYODIDE_URL}pyodide.js`; script.onload = resolve; script.onerror = reject; document.head.append(script); });
    state.pyodide = await window.loadPyodide({ indexURL: PYODIDE_URL });
    await state.pyodide.runPythonAsync(parserCode);
    elements.status.textContent = "Local parser ready";
  } catch (error) {
    elements.status.textContent = "Parser unavailable";
    console.error(error);
  }
}
async function importFiles(files) {
  if (!state.pyodide) { elements.summary.textContent = "The local parser is still loading; please try again in a moment."; return; }
  const added = [];
  for (const file of files) {
    try {
      elements.status.textContent = `Parsing ${file.name}…`;
      state.pyodide.globals.set("gb_text", await file.text());
      const decoded = JSON.parse(await state.pyodide.runPythonAsync("parse_genbank(gb_text)"));
      if (!decoded.length) throw new Error("No GenBank FEATURES section was found");
      decoded.forEach((genome, i) => added.push({ ...genome, name: decoded.length > 1 ? `${genome.name} (${i + 1})` : genome.name, fileName: file.name }));
    } catch (error) { console.warn(error); elements.summary.textContent = `${file.name}: ${error.message}. Other valid files were retained.`; }
  }
  state.genomes.push(...added);
  elements.summary.textContent = state.genomes.length ? `${state.genomes.length} genome${state.genomes.length === 1 ? "" : "s"} loaded from ${new Set(state.genomes.map(g => g.fileName)).size} file(s).` : "No valid GenBank genomes loaded.";
  elements.status.textContent = "Local parser ready";
  render();
}

elements.files.addEventListener("change", event => importFiles([...event.target.files]));
elements.render.addEventListener("click", render);
[...document.querySelectorAll('input[name="view"]'), elements.theme, elements.start, elements.labels, elements.strands].forEach(input => input.addEventListener("change", render));
elements.reset.addEventListener("click", () => { state.genomes = []; elements.files.value = ""; elements.summary.textContent = "No genomes loaded."; elements.stats.innerHTML = ""; elements.plot.className = "plot empty"; elements.plot.innerHTML = "<div><span class=\"empty-icon\">◎</span><p>Your mitochondrial maps will appear here.</p></div>"; elements.title.textContent = "Waiting for GenBank files"; elements.note.textContent = "Load one or more annotated genomes to begin."; elements.download.disabled = true; });
elements.download.addEventListener("click", () => { const svg = elements.plot.querySelector("svg"); if (!svg) return; const blob = new Blob([new XMLSerializer().serializeToString(svg)], { type:"image/svg+xml;charset=utf-8" }); const link = Object.assign(document.createElement("a"), { href:URL.createObjectURL(blob), download:`pyvam-${viewMode()}-maps.svg` }); link.click(); URL.revokeObjectURL(link.href); });

loadPyodide();
