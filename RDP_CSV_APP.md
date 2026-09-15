# RDP CSV → Pretty Pictures

This repo now includes a focused local Dash app for turning RDP-exported recombination event tables into interactive figures.

## Run it locally

From the repository directory:

```bash
conda create -n rdp-dash python=3.11 -y
conda activate rdp-dash
pip install -r requirements.txt
python rdp_csv_visualizer.py
```

Then open:

**http://127.0.0.1:8050**

If port 8050 is already being used by the main `app.py`, stop that process first or change the port at the bottom of `rdp_csv_visualizer.py`.

## Input

Upload the CSV/TSV results exported by RDP after a recombination scan. The parser recognizes common variants of:

- Recombinant / recombinant sequence
- Major parent
- Minor parent
- Beginning/start breakpoint
- Ending/end breakpoint
- Event number
- Per-method result or p-value columns such as RDP, GENECONV, BootScan, MaxChi, Chimaera, SiScan, and 3Seq

Column names are normalized, so spaces, punctuation, capitalization, and several common aliases are accepted.

## Figures

### Event map
Each recombinant is shown as a genome/alignment track. Detected recombinant regions are drawn as colored blocks positioned at their reported breakpoints and colored by minor parent.

### Parent flow
A Sankey-style diagram summarizes major-parent → minor-parent → recombinant relationships across events.

### Breakpoint density
Start and end breakpoints are pooled into a genome-wide density histogram to make recurrent breakpoint regions visually obvious.

### Method support
A heatmap shows which RDP methods support each detected event when recognizable per-method fields are present in the export.

### Clean table
A sortable/filterable normalized event table that can also be exported as CSV.

## Exporting figures

Hover over a Plotly figure and click the camera icon to export a 2× PNG suitable for slides and teaching material.

> This app visualizes results reported by RDP. It does not rerun recombination detection or independently validate an event.
