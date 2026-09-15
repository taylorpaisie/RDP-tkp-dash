# RDP5 Folder Visualizer

This mode is designed for RDP5 workflows where project files are written to a local directory and should appear in the browser automatically.

## What it does

- Watches a local folder for `*.rdp5` files.
- Rescans every 3 seconds by default.
- Automatically selects the newest project when `Automatically open newest project` is enabled.
- Reads the selected RDP5 project directly from disk; there is no browser upload step.
- Extracts the sequence labels and embedded alignment structures currently supported by `app.py`.
- Renders alignment divergence, genome-wide signal profiles, subtype composition, country composition, and a searchable label table.

The supplied VEME example project was used to validate the current parser. It contains 274 sequence labels and 274 embedded sequence runs of approximately 9.56 kb.

## Run

```bash
pip install -r requirements.txt
python rdp5_folder_app.py
```

Open:

```text
http://127.0.0.1:8050
```

The app intentionally binds to `127.0.0.1` because it can read a configured local folder.

## Default folder

If no environment variable is supplied, the app creates and watches:

```text
./rdp5_projects
```

Drop `.rdp5` files there and they will appear automatically.

## Watch an existing RDP output directory

### Linux / macOS / WSL

```bash
export RDP5_DIR="/path/to/RDP/output"
python rdp5_folder_app.py
```

### PowerShell

```powershell
$env:RDP5_DIR = "C:\path\to\RDP\output"
python rdp5_folder_app.py
```

You can also type another folder path into the web page and click **Use folder**.

## Poll interval

The default is 3000 ms. Override it with:

```bash
export RDP5_POLL_MS=1000
```

Avoid very aggressive polling when RDP is writing large project files.

## Current RDP5 decoding boundary

The RDP5 binary format is not publicly documented. The example file lets us reliably validate:

- the RDP5 signature,
- sequence-label records,
- embedded aligned nucleotide sequences.

The app does not yet claim to decode RDP5's persisted recombination-event structures directly from the binary file. Those should only be added when their binary layout has been validated against known RDP5 projects/results.
