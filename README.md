# RDP-tkp Visualizer

A Dash app for exploring RDP teaching files and exported recombination-event results. It complements **[nextRDP Web](https://murrellgroup.github.io/nextRDPweb/)**; it does not rerun or replace RDP.

## What it visualizes

### FASTA alignments

- sequence count and aligned length;
- per-sequence gap and ambiguity rates;
- sliding-window variable-site counts;
- sliding-window nucleotide entropy; and
- sortable sequence-level QC.

### Exported RDP results

Upload a CSV or TSV containing at least recombinant, start, and end columns. Common column-name variants are recognized. The app shows:

- candidate breakpoint intervals;
- events by recombinant and method;
- parental labels and p-values when present; and
- a sortable/filterable event table.

### RDP5 project files

The app validates the `RDP5 Project File` signature, extracts the embedded alignment and fixed-width sequence-label inventory, and renders an RDP-style alignment overview plus a genome-wide divergence/gap profile. Label metadata remain available in the sortable table. The binary `.rdp5` event format is not publicly documented, so the app deliberately does **not** claim event-level parsing from project binaries. Export results to CSV/TSV from RDP for full event visualization.

The **Similarity scan** tab provides a BootScan-style teaching view. Select a query or suspected recombinant, compare it with up to four candidate parents, and adjust the sliding-window and step sizes. Curves show pairwise nucleotide identity at sites called A/C/G/T in both sequences. Similarity switches identify regions worth investigating; they are not independently interpreted as recombination calls.

Below the curves, a **candidate recombination mosaic** colors each window by its closest selected parent and marks parent-switch transitions. Hover text reports the winning parent and its percentage-point lead over the next-closest parent. These transitions are exploratory candidates, not confirmed RDP breakpoints.

## VEME lesson files tested

| File | Type | Observed structure |
|---|---|---|
| Exercise 1 alignment.fas | FASTA alignment | 25 sequences × 9,594 sites |
| Exercise 2 alignment.fas | FASTA alignment | 50 sequences × 8,108 sites |
| Exercise 3 RDP project file.rdp5 | RDP5 project | 34.8 MB binary project with a valid RDP5 signature |

The course files are not committed here because the repository is public and the large project file may contain bundled reference data. Students can upload the copies distributed during the lesson.

## Install

```bash
git clone https://github.com/taylorpaisie/RDP-tkp-dash.git
cd RDP-tkp-dash
conda create -n rdp-tkp python=3.11 pip -y
conda activate rdp-tkp
pip install -r requirements.txt
python app.py
```

Open <http://127.0.0.1:8050>.

## Docker

```bash
docker build -t rdp-tkp .
docker run --rm -p 8050:8050 rdp-tkp
```

## Recommended lesson flow

1. Upload Exercise 1 and discuss alignment QC and localized variability.
2. Upload Exercise 2 and compare diversity/gap profiles with Exercise 1.
3. Upload the Exercise 3 project to inspect its project/sequence inventory.
4. Open [nextRDP Web](https://murrellgroup.github.io/nextRDPweb/) for recombination detection and evidence review.
5. Export an event table and return to this app for breakpoint and method-level visualization.

## Upstream resources

- [nextRDP Web](https://murrellgroup.github.io/nextRDPweb/)
- [nextRDP Web source](https://github.com/MurrellGroup/nextRDPweb)
- [nextRDP core](https://github.com/MurrellGroup/nextRDP-core)
- [RDP home page](https://web.cbio.uct.ac.za/~darren/rdp.html)
- [RDP5 paper](https://doi.org/10.1093/ve/veaa087)

## Disclaimer

Teaching and exploratory visualization only. Candidate events require alignment review, method-specific evidence assessment, phylogenetic context, and biological interpretation.
