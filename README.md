# RDP-tkp — VEME 2026 practical

Hands-on materials for the VEME recombination-detection lesson led by Darren Martin, with Taylor Paisie assisting.

## Start here

1. Open the [student lesson](LESSON.md).
2. Download the [synthetic teaching alignment](data/veme_synthetic_recombination.fasta).
3. Launch **[nextRDP Web](https://murrellgroup.github.io/nextRDPweb/)**.
4. Record observations in the [student worksheet](WORKSHEET.md).

Instructors: use the [run sheet](INSTRUCTOR_GUIDE.md) and [answer key](ANSWER_KEY.md).

## Learning objectives

By the end of the practical, participants should be able to:

- explain why recombination can mislead phylogenetic inference;
- recognize a recombinant sequence and candidate parental lineages;
- compare support across multiple detection methods;
- inspect breakpoint evidence and regional phylogenies;
- distinguish an automated candidate event from a reviewed biological conclusion; and
- export a reproducible result record.

## Why WebRDP

nextRDP Web is a browser interface for the source-faithful, WebAssembly-compatible nextRDP core. It provides RDP, GENECONV, MaxChi, CHIMAERA, 3SEQ, BootScan, and SISCAN discovery lanes, evidence plots, breakpoint alignments, regional trees, PHYLPRO profiles, review controls, and exports. Alignment data remain inside the browser.

> The upstream developers currently label the BootScan and SISCAN lanes as source-shaped but unvalidated pending full cyclic integration. Treat them as supporting/exploratory evidence during this lesson.

## Repository map

- `LESSON.md` — student-facing practical
- `WORKSHEET.md` — observations and interpretation prompts
- `INSTRUCTOR_GUIDE.md` — timing, teaching cues, troubleshooting, and debrief
- `ANSWER_KEY.md` — expected qualitative findings for the synthetic example
- `data/veme_synthetic_recombination.fasta` — aligned, synthetic teaching data

## Upstream resources

- [nextRDP Web](https://murrellgroup.github.io/nextRDPweb/)
- [nextRDP Web source](https://github.com/MurrellGroup/nextRDPweb)
- [nextRDP core](https://github.com/MurrellGroup/nextRDP-core)
- [RDP home page](https://web.cbio.uct.ac.za/~darren/rdp.html)
- [RDP5 paper](https://doi.org/10.1093/ve/veaa087)

## Citation

Martin DP, Varsani A, Roumagnac P, Botha G, Maslamoney S, Schwab T, Kelz Z, Kumar V, and Murrell B. (2021). RDP5: a computer program for analyzing recombination in, and removing signals of recombination from, nucleotide sequence datasets. *Virus Evolution*, 7, veaa087. https://doi.org/10.1093/ve/veaa087

## Data note

The included alignment is entirely synthetic and designed only for teaching. It contains no patient, outbreak, or unpublished research data.
