from __future__ import annotations

import base64
import io
import subprocess
import tempfile
from pathlib import Path

import dash_bootstrap_components as dbc
import pandas as pd
import plotly.express as px
from dash import Dash, Input, Output, State, callback, dash_table, dcc, html, no_update


METHODS = ["geneconv", "bootscan", "maxchi", "siscan", "chimaera", "threeseq", "rdp"]
COLUMNS = ["Method", "Start", "End", "Recombinant", "Parent1", "Parent2", "Pvalue"]

DEMO = pd.DataFrame(
    [
        ["RDP", 6, 504, "X64860", "X64866", "X64869", 6.87e-7],
        ["3Seq", 202, 787, "X64869", "X64860", "X64866", 5.98e-10],
        ["MaxChi", 439, 482, "X64860", "X64866", "X64869", 4.04e-2],
        ["Chimaera", 179, 265, "X64866", "X64869", "X64873", 4.70e-3],
    ],
    columns=COLUMNS,
)


def decode_upload(contents: str) -> str:
    _, encoded = contents.split(",", 1)
    return base64.b64decode(encoded).decode("utf-8")


def validate_fasta(text: str) -> tuple[int, int]:
    names, sequences, current = [], [], []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current:
                sequences.append("".join(current))
                current = []
            names.append(line[1:].strip())
        elif not names:
            raise ValueError("Sequence data appears before the first FASTA header.")
        else:
            current.append(line)
    if current:
        sequences.append("".join(current))
    if len(names) < 3 or len(names) != len(sequences):
        raise ValueError("Provide at least three complete FASTA records.")
    lengths = {len(s) for s in sequences}
    if len(lengths) != 1:
        raise ValueError("Sequences must already be aligned and have equal lengths.")
    return len(names), lengths.pop()


def run_openrdp(fasta: str, methods: list[str]) -> pd.DataFrame:
    with tempfile.TemporaryDirectory(prefix="openrdp_dash_") as tmp:
        infile, outfile = Path(tmp) / "alignment.fasta", Path(tmp) / "results.csv"
        infile.write_text(fasta, encoding="utf-8")
        cmd = ["openrdp", str(infile), "-o", str(outfile), "-q"]
        if methods:
            cmd.extend(["-m", *methods])
        try:
            completed = subprocess.run(cmd, capture_output=True, text=True, timeout=1800, check=False)
        except FileNotFoundError as exc:
            raise RuntimeError("OpenRDP is not installed. Run: pip install -r requirements.txt") from exc
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError("The scan exceeded the 30-minute safety limit.") from exc
        if completed.returncode != 0:
            detail = (completed.stderr or completed.stdout or "Unknown OpenRDP error").strip()
            raise RuntimeError(detail[-1200:])
        if not outfile.exists():
            raise RuntimeError("OpenRDP completed without creating its CSV output.")
        frame = pd.read_csv(outfile)
        frame.columns = [str(c).strip() for c in frame.columns]
        return frame


app = Dash(__name__, external_stylesheets=[dbc.themes.FLATLY], title="OpenRDP Dash")
server = app.server

app.layout = dbc.Container(
    [
        dcc.Store(id="results-store"),
        dbc.Row(
            dbc.Col(
                [html.H1("OpenRDP Dash", className="display-5 fw-bold mb-1"),
                 html.P("Upload an aligned FASTA, run recombination scans, and inspect candidate breakpoint intervals.", className="lead text-secondary")],
                className="py-4",
            )
        ),
        dbc.Row(
            [
                dbc.Col(
                    dbc.Card(dbc.CardBody([
                        html.H4("1 · Analysis", className="card-title"),
                        dcc.Upload(id="upload", children=html.Div(["Drop an aligned FASTA here or ", html.A("browse")]),
                                   className="border rounded p-4 text-center mb-3", multiple=False),
                        html.Div(id="file-status", className="small mb-3"),
                        dbc.Label("Detection methods"),
                        dbc.Checklist(id="methods", options=[{"label": m, "value": m} for m in METHODS], value=METHODS,
                                      inline=True, className="mb-3"),
                        dbc.Button("Run OpenRDP", id="run", color="primary", className="me-2"),
                        dbc.Button("Load demo", id="demo", color="secondary", outline=True),
                        html.Div(id="run-status", className="mt-3"),
                    ])), md=4,
                ),
                dbc.Col(
                    dbc.Card(dbc.CardBody([
                        html.H4("2 · Results", className="card-title"),
                        dbc.Row([
                            dbc.Col([dbc.Label("Methods"), dcc.Dropdown(id="method-filter", multi=True)], md=7),
                            dbc.Col([dbc.Label("Maximum p-value"), dbc.Input(id="p-filter", type="number", value=0.05, min=0, max=1, step=0.01)], md=5),
                        ], className="mb-3"),
                        dcc.Graph(id="breakpoint-chart", config={"displaylogo": False}),
                        dash_table.DataTable(id="results-table", page_size=10, sort_action="native", filter_action="native",
                                             style_table={"overflowX": "auto"}, style_cell={"fontFamily": "system-ui", "fontSize": 13, "padding": "8px"}),
                        dbc.Button("Download filtered CSV", id="download-button", color="success", outline=True, className="mt-3"),
                        dcc.Download(id="download"),
                    ])), md=8,
                ),
            ], className="g-3 pb-4"
        ),
        html.P(["Experimental interface powered by ", html.A("OpenRDP", href="https://github.com/PoonLab/OpenRDP", target="_blank"),
                ". Review upstream licenses before commercial use."], className="small text-secondary text-center"),
    ], fluid="xl"
)


@callback(Output("file-status", "children"), Input("upload", "contents"), State("upload", "filename"))
def inspect_file(contents, filename):
    if not contents:
        return "No alignment selected."
    try:
        count, length = validate_fasta(decode_upload(contents))
        return dbc.Alert(f"{filename}: {count} sequences × {length:,} sites", color="success", className="py-2")
    except Exception as exc:
        return dbc.Alert(str(exc), color="danger", className="py-2")


@callback(
    Output("results-store", "data"), Output("run-status", "children"),
    Input("run", "n_clicks"), Input("demo", "n_clicks"),
    State("upload", "contents"), State("methods", "value"), prevent_initial_call=True,
)
def execute(run_clicks, demo_clicks, contents, methods):
    from dash import ctx
    if ctx.triggered_id == "demo":
        return DEMO.to_dict("records"), dbc.Alert("Demo results loaded.", color="info", className="py-2")
    if not contents:
        return no_update, dbc.Alert("Upload an aligned FASTA first.", color="warning", className="py-2")
    try:
        fasta = decode_upload(contents)
        validate_fasta(fasta)
        frame = run_openrdp(fasta, methods or [])
        return frame.to_dict("records"), dbc.Alert(f"Scan complete: {len(frame):,} events.", color="success", className="py-2")
    except Exception as exc:
        return no_update, dbc.Alert(str(exc), color="danger", className="py-2")


def filtered_frame(records, selected, pmax):
    frame = pd.DataFrame(records or [])
    if frame.empty:
        return frame
    pcol = next((c for c in frame.columns if c.lower().replace(" ", "") in {"pvalue", "p-value"}), None)
    mcol = next((c for c in frame.columns if c.lower() == "method"), None)
    if selected and mcol:
        frame = frame[frame[mcol].astype(str).isin(selected)]
    if pmax is not None and pcol:
        frame = frame[pd.to_numeric(frame[pcol], errors="coerce") <= float(pmax)]
    return frame


@callback(
    Output("method-filter", "options"), Output("method-filter", "value"),
    Output("breakpoint-chart", "figure"), Output("results-table", "data"), Output("results-table", "columns"),
    Input("results-store", "data"), Input("method-filter", "value"), Input("p-filter", "value"),
)
def render(records, selected, pmax):
    base = pd.DataFrame(records or [])
    if base.empty:
        fig = px.scatter(title="Run a scan or load the demo to see results")
        fig.update_layout(template="plotly_white")
        return [], [], fig, [], []
    mcol = next((c for c in base.columns if c.lower() == "method"), "Method")
    methods = sorted(base[mcol].astype(str).unique())
    options = [{"label": m, "value": m} for m in methods]
    active = selected if selected is not None and set(selected).issubset(methods) else methods
    frame = filtered_frame(records, active, pmax)
    start = next((c for c in frame.columns if c.lower() == "start"), "Start")
    end = next((c for c in frame.columns if c.lower() == "end"), "End")
    recombinant = next((c for c in frame.columns if c.lower() == "recombinant"), "Recombinant")
    if frame.empty:
        fig = px.scatter(title="No events match the current filters")
    else:
        chart = frame.copy()
        chart["Length"] = pd.to_numeric(chart[end]) - pd.to_numeric(chart[start]) + 1
        fig = px.bar(chart, x="Length", y=recombinant, base=start, color=mcol, orientation="h",
                     hover_data=[start, end], title="Candidate recombination intervals")
    fig.update_layout(template="plotly_white", xaxis_title="Alignment position", yaxis_title="Recombinant")
    return options, active, fig, frame.to_dict("records"), [{"name": c, "id": c} for c in frame.columns]


@callback(Output("download", "data"), Input("download-button", "n_clicks"),
          State("results-store", "data"), State("method-filter", "value"), State("p-filter", "value"), prevent_initial_call=True)
def download_results(_, records, selected, pmax):
    return dcc.send_data_frame(filtered_frame(records, selected, pmax).to_csv, "openrdp_results.csv", index=False)


if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=8050)

