/* Cancellable PyVAM figure renderer. Heavy Matplotlib work stays off the UI thread. */
const PYODIDE_URL = "https://cdn.jsdelivr.net/pyodide/v0.27.5/full/";
let pyodide = null;
let initialized = false;
const respond = (id, result) => postMessage({ id, ok: true, result });
const fail = (id, error) => postMessage({ id, ok: false, error: String(error?.message || error) });
async function initialize(runtime) {
  if (initialized) return;
  importScripts(PYODIDE_URL + "pyodide.js");
  pyodide = await loadPyodide({ indexURL: PYODIDE_URL });
  await pyodide.loadPackage(["matplotlib", "biopython"]);
  pyodide.FS.mkdirTree("/home/pyvam");
  for (const name of ["config.py", "parserGB.py", "drawMT.py"]) {
    const response = await fetch("./pyvam/" + name);
    if (!response.ok) throw new Error("PyVAM source unavailable: " + name);
    pyodide.FS.writeFile("/home/pyvam/" + name, await response.text());
  }
  pyodide.FS.writeFile("/home/pyvam/__init__.py", "from .parserGB import get_features\nfrom .drawMT import draw_circos_MT, draw_linear_MT, draw_linear_MT_nonproportional\n");
  await pyodide.runPythonAsync("import sys; sys.path.insert(0, '/home')");
  await pyodide.runPythonAsync(runtime);
  initialized = true;
}
function ensureParent(path) {
  const components = path.split("/").filter(Boolean);
  let current = "";
  for (const component of components.slice(0, -1)) {
    current += "/" + component;
    try { pyodide.FS.mkdir(current); } catch (error) { if (error?.code !== "EEXIST") throw error; }
  }
}
self.onmessage = async ({ data }) => {
  const { id, type, payload } = data;
  try {
    if (type === "init") {
      await initialize(payload.runtime);
      respond(id, true);
    } else if (type === "render") {
      if (!initialized) throw new Error("Renderer not initialized");
      for (const file of payload.files) {
        ensureParent(file.path);
        pyodide.FS.writeFile(file.path, file.content);
      }
      for (const [key, value] of Object.entries(payload.globals)) pyodide.globals.set(key, value);
      const result = await pyodide.runPythonAsync("pyvam_render(web_paths_json, web_view, web_theme, web_start, web_labels, web_custom_colors, web_options)");
      respond(id, result);
    } else if (type === "export") {
      if (!initialized) throw new Error("Renderer not initialized");
      const result = await pyodide.runPythonAsync("pyvam_export(" + JSON.stringify(payload.format) + ")");
      respond(id, result);
    } else {
      throw new Error("Unsupported worker command: " + type);
    }
  } catch (error) { fail(id, error); }
};
