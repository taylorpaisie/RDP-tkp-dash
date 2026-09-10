# RDP-tkp

A compact workflow companion for **[nextRDP Web](https://murrellgroup.github.io/nextRDPweb/)**, the browser-based RDP implementation recommended by Darren Martin.

## Launch WebRDP

**[Open nextRDP Web →](https://murrellgroup.github.io/nextRDPweb/)**

No Python, Miniconda, Dash server, or local OpenRDP installation is required. Analysis runs in the browser through a WebAssembly-compatible, source-faithful RDP core, and alignment data stay in the browser.

## What WebRDP currently provides

- RDP, GENECONV, MaxChi, CHIMAERA, 3SEQ, BootScan, and SISCAN discovery lanes
- Method-specific evidence plots
- Breakpoint alignment review
- Regional trees and PHYLPRO profiles
- Review controls and export views

> BootScan and SISCAN are currently described by the developers as source-shaped but unvalidated until full cyclic integration is completed. Treat those lanes accordingly.

## Suggested workflow

1. Prepare a nucleotide multiple-sequence alignment.
2. Open WebRDP and load the alignment.
3. Choose the appropriate discovery methods and settings.
4. Review supported events across methods rather than accepting raw calls automatically.
5. Inspect breakpoint alignments, evidence plots, and regional phylogenies.
6. Export the reviewed results and record the WebRDP/core version used.

## Upstream projects

- [nextRDP Web application](https://github.com/MurrellGroup/nextRDPweb)
- [nextRDP core](https://github.com/MurrellGroup/nextRDP-core)
- [RDP home page](https://web.cbio.uct.ac.za/~darren/rdp.html)
- [RDP5 paper](https://doi.org/10.1093/ve/veaa087)

## Citation

Martin DP, Varsani A, Roumagnac P, Botha G, Maslamoney S, Schwab T, Kelz Z, Kumar V, and Murrell B. (2021). RDP5: a computer program for analyzing recombination in, and removing signals of recombination from, nucleotide sequence datasets. *Virus Evolution*, 7, veaa087. https://doi.org/10.1093/ve/veaa087

## Repository history

The initial commit explored a separate Python Dash wrapper around OpenRDP. That approach was retired after consultation with Darren Martin in favor of the actively developed, source-faithful WebRDP application. The original prototype remains available in Git history.
