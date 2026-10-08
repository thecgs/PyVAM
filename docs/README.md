# PyVAM Web

This static web companion loads Pyodide from its CDN, then loads bundled PyVAM source and calls PyVAM's own `draw_circos_MT`, `draw_linear_MT`, and `draw_linear_MT_nonproportional` functions. Selected GenBank files stay in the visitor's browser; they are written only to Pyodide's in-memory virtual filesystem.

It supports selecting multiple files, circular maps, proportional linear maps, gene-order comparison, config-defined themes, config-defined alignment markers, gene-level colour overrides, and SVG export. The image is a Matplotlib SVG emitted by PyVAM, rather than a JavaScript reimplementation.

## Deploy with GitHub Pages

The included workflow deploys this directory after a push to `main`. In the repository's **Settings → Pages**, set **Source** to **GitHub Actions** once. For a fork using another default branch, update `.github/workflows/pages.yml` accordingly.
