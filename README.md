# RDP-tkp Visualizer

An interactive Dash app for teaching, exploring, and presenting recombination evidence from RDP lesson files. It complements [RDP](https://web.cbio.uct.ac.za/~darren/rdp.html) and [nextRDP Web](https://murrellgroup.github.io/nextRDPweb/); it does not rerun or replace their statistical detection methods.

## Features

### RDP5 project overview

Upload an `.rdp5` project to extract its embedded alignment and sequence labels. The overview provides:

- an RDP-style alignment heatmap showing divergence from the consensus;
- a genome-wide mean-divergence profile;
- gap and uncalled-site coverage;
- sequence, site, subtype, and country-code summaries; and
- a searchable sequence-metadata table.

The heatmap divides the alignment into 120 windows. For each sequence and window, divergence is calculated as:

```text
100 × called nucleotides different from the consensus / called nucleotides
```

Only A, C, G, and T contribute to the consensus and denominator. Large projects display up to 72 representative tracks spanning the observed divergence range.

### Similarity scan

The **Similarity scan** tab provides a BootScan-style teaching view:

1. Choose a query or suspected recombinant sequence.
2. Select up to four candidate parents.
3. Set a sliding window from 100–800 nucleotides.
4. Choose a step of 25, 50, 100, or 200 nucleotides.

The plot shows pairwise nucleotide identity between the query and each selected parent. Sites are included only when both sequences contain A, C, G, or T. The identity axis automatically zooms to the informative range.

Candidate parents keep the same color in every panel using a colorblind-friendly palette:

| Parent | Color |
|---|---|
| 1 | Blue |
| 2 | Orange |
| 3 | Teal |
| 4 | Magenta |

### Candidate recombination mosaic

The mosaic summarizes similarity switching along the genome:

- each segment is colored by the closest selected parent;
- adjacent windows assigned to the same parent are merged;
- vertical lines indicate closest-parent transitions;
- a five-window categorical smoother suppresses isolated parent flicker; and
- the confidence panel shows the identity lead over the second-closest parent.

Confidence is expressed in percentage points. A value near zero indicates that the two best parents are effectively tied. A larger, sustained lead makes the displayed parental assignment more distinct, but it does not establish statistical support for recombination.

Hover over a segment to inspect its coordinates, closest parent, mean identity, and mean confidence.

### FASTA alignment QC

FASTA inputs provide:

- sequence count and aligned length;
- per-sequence gap and ambiguity rates;
- sliding-window variable-site counts;
- sliding-window nucleotide entropy; and
- sortable sequence-level QC.

### Exported RDP event tables

CSV and TSV tables must contain recombinant, start, and end fields. Common column-name variants are accepted. When present, method, parental-label, and p-value fields are also retained. The app displays:

- candidate breakpoint intervals;
- events by recombinant and method;
- parental labels and p-values; and
- a sortable and filterable event table.

## Important interpretation note

The `.rdp5` project format is a proprietary binary format. This app extracts the embedded alignment and sequence names, but it does **not** claim to decode RDP's stored event calls, corrected p-values, or method support from the binary file.

The similarity scan and mosaic are exploratory alignment visualizations. A change in the closest parent is a region to investigate—not a confirmed breakpoint. Recombination conclusions should incorporate RDP results, multiple-method support, breakpoint uncertainty, phylogenetic evidence, and biological context.

For event-level visualizations, export a CSV or TSV results table from RDP and upload it separately.

## VEME 2026 lesson files tested

| Exercise | Input | Validated structure |
|---|---|---:|
| Exercise 1 | FASTA and RDP5 | 25 sequences × 9,594 sites |
| Exercise 2 | FASTA and RDP5 | 50 sequences × 8,108 sites |
| Exercise 3 | RDP5 | 274 sequences × 9,556 sites; 34.8 MB project |

All three RDP5 projects were tested through label extraction, alignment parsing, overview rendering, similarity scanning, mosaic generation, and metadata-table rendering.

The course files are not committed because the repository is public and the projects may bundle reference data. Students should upload the files distributed with the lesson.

## Installation

### Conda

```bash
git clone https://github.com/taylorpaisie/RDP-tkp-dash.git
cd RDP-tkp-dash
conda create -n rdp-tkp python=3.11 pip -y
conda activate rdp-tkp
pip install -r requirements.txt
python app.py
```

Open <http://127.0.0.1:8050>.

If `conda` is unavailable inside WSL, install Miniconda in WSL first. A Conda installation on Windows is not automatically available inside the Linux environment.

### Docker

```bash
docker build -t rdp-tkp .
docker run --rm -p 8050:8050 rdp-tkp
```

## Suggested classroom workflow

1. Upload Exercise 1 and review alignment quality, variable regions, and gaps.
2. Open the similarity scan and choose a suspected recombinant and candidate parents.
3. Change the window size to demonstrate the resolution-versus-noise tradeoff.
4. Interpret the mosaic alongside its confidence panel, emphasizing weak near-ties.
5. Repeat with Exercise 2 and compare the genome-wide signal.
6. Use Exercise 3 to demonstrate scaling to hundreds of sequences.
7. Review the formal event evidence in RDP or nextRDP Web.
8. Upload an exported event table for breakpoint and method-level summaries.

## Upstream resources

- [RDP home page](https://web.cbio.uct.ac.za/~darren/rdp.html)
- [RDP5 paper](https://doi.org/10.1093/ve/veaa087)
- [nextRDP Web](https://murrellgroup.github.io/nextRDPweb/)
- [nextRDP Web source](https://github.com/MurrellGroup/nextRDPweb)
- [nextRDP core](https://github.com/MurrellGroup/nextRDP-core)

## Disclaimer

For teaching and exploratory visualization only. Inspect primary evidence before making biological conclusions.
