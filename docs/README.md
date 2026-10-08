# PyVAM Web

This is a no-build, static web companion for PyVAM. It loads Pyodide from its CDN and parses selected GenBank files in the visitor's browser; sequence files are not uploaded by this application.

It supports selecting multiple files (and multi-record GenBank files), circular maps, linear maps, gene-order comparison, themes, rotation at a chosen marker, and SVG export.

## Deploy with GitHub Pages

The included workflow deploys this directory after a push to `main`. In the repository's **Settings → Pages**, set **Source** to **GitHub Actions** once. For a fork using another default branch, update `.github/workflows/pages.yml` accordingly.
