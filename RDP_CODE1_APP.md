# RDP CSV Code 1 local viewer

`rdp_code1_app.py` is a focused local Dash app for the **RDP CSV Code 1** event-plot export.

It is designed for the workflow where RDP writes CSV exports to a directory and the browser app should notice new files automatically, rather than requiring the user to upload each file manually.

## What it reads

The parser is based on the supplied RDP v5.93 Code 1 example and reads:

- gene start/end coordinates;
- reading frame and orientation;
- event number and axis labels;
- beginning and ending breakpoint sites;
- 95% and 99% breakpoint confidence intervals;
- RDP plot colors and pairwise-comparison roles; and
- the position-by-position plot data.

## What it renders

The browser view contains:

- a gene-map track;
- the three RDP pairwise-identity curves;
- shaded 95% and 99% breakpoint confidence intervals;
- vertical breakpoint-site markers;
- event, breakpoint, gene-count, and plotted-position summary cards; and
- Plotly PNG export from the figure toolbar.

The pairwise values in the supplied export are integer-scaled. The viewer divides all three curves by the largest observed value in the file so the displayed y-axis matches RDP's 0–1 pairwise-identity presentation.

## Run locally

```bash
conda create -n rdp-tkp python=3.11 pip -y
conda activate rdp-tkp
pip install -r requirements.txt
python rdp_code1_app.py
```

Then open:

```text
http://127.0.0.1:8050
```

## Folder watching

By default the app creates and watches:

```text
./rdp_code1_csv/
```

Place RDP Code 1 CSV exports in that directory. The app scans every 3 seconds and, by default, switches to the newest CSV automatically.

To point it at RDP's existing export directory, set `RDP_CODE1_DIR` before starting the app.

### PowerShell

```powershell
$env:RDP_CODE1_DIR = "C:\path\to\RDP\exports"
python rdp_code1_app.py
```

### Bash / WSL

```bash
export RDP_CODE1_DIR="/path/to/RDP/exports"
python rdp_code1_app.py
```

The polling interval can also be changed. For example, to scan every 5 seconds:

```powershell
$env:RDP_CODE1_POLL_MS = "5000"
```

## Scope

This first parser intentionally targets **CSV Code 1**, because that is the concrete RDP export supplied for the desired event plot. Other RDP CSV codes should be implemented as explicit parsers once example files are available rather than guessed from undocumented formats.

The web server binds to `127.0.0.1` so a process that can read local directories is not exposed to the local network by default.
