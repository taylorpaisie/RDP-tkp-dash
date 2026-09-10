# OpenRDP Dash

An experimental web interface for running and exploring recombination scans with [OpenRDP](https://github.com/PoonLab/OpenRDP). It is inspired by the broader analysis workflow in [RDP5](https://web.cbio.uct.ac.za/~darren/rdp.html), but it is not affiliated with or a replacement for RDP5.

## Features

- Upload aligned nucleotide sequences in FASTA format
- Validate names and alignment lengths before analysis
- Select any of OpenRDP's seven methods
- Run OpenRDP through its command-line interface
- Filter results by method and p-value
- Explore breakpoint intervals in an interactive Plotly chart
- Download filtered results as CSV
- Load a built-in demo result without installing OpenRDP

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python app.py
```

Open <http://127.0.0.1:8050>.

OpenRDP is installed directly from its upstream GitHub repository. Some bundled third-party methods have non-commercial or academic-use license restrictions; review the upstream licenses before use. Uploaded files are written only to a temporary directory and removed after each run.

## Docker

```bash
docker build -t openrdp-dash .
docker run --rm -p 8050:8050 openrdp-dash
```

## Scope

This first version is a thin, reproducible web layer over the OpenRDP CLI. A useful next phase would expose method-specific configuration, add alignment previews, persist run manifests, and test output compatibility against curated RDP5 examples.

## Disclaimer

Research software only. Results should be reviewed with appropriate controls and biological context.

