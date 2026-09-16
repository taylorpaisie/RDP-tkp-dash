from __future__ import annotations

import os
import re
from pathlib import Path

import dash_bootstrap_components as dbc
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from dash import Dash, Input, Output, State, callback, dcc, html, no_update

POLL_MS = int(os.getenv("RDP_CODE1_POLL_MS", "3000"))
DEFAULT_DIR = Path(os.getenv("RDP_CODE1_DIR", Path(__file__).with_name("rdp_code1_csv"))).expanduser()
if "RDP_CODE1_DIR" not in os.environ:
    DEFAULT_DIR.mkdir(parents=True, exist_ok=True)

COLORS = {"yellow": "#e2d51d", "purple": "#b948bd", "green": "#20a56a", "blue": "#2ea8c7"}
FALLBACK = ["#e2d51d", "#b948bd", "#20a56a"]
FRAME_Y = {1: 2.5, 4: 2.5, 2: 1.5, 5: 1.5, 3: 0.5, 6: 0.5}


def folder_path(value: str) -> Path:
    p = Path(value).expanduser().resolve()
    if not p.is_dir():
        raise ValueError(f"Folder does not exist: {p}")
    return p


def discover(folder: str) -> list[Path]:
    return sorted(folder_path(folder).glob("*.csv"), key=lambda p: p.stat().st_mtime_ns, reverse=True)


def parse_pair(values: list[str]):
    try:
        return int(values[0]), int(values[1])
    except (ValueError, IndexError):
        return None, None


def parse_code1(path: Path) -> dict:
    lines = path.read_text(encoding="utf-8-sig", errors="replace").splitlines()
    if not lines or not lines[0].startswith("Gene start"):
        raise ValueError("Not an RDP CSV Code 1 export.")

    try:
        split = next(i for i, line in enumerate(lines) if not line.strip())
        plot = lines.index("Plot data")
    except (StopIteration, ValueError) as exc:
        raise ValueError("RDP Code 1 sections were not found.") from exc

    genes = []
    for line in lines[1:split]:
        parts = [x.strip() for x in line.split(",")]
        if len(parts) < 4:
            continue
        try:
            start, end, frame, orient = map(int, parts[:4])
            genes.append({"start": start, "end": end, "frame": frame, "orientation": orient})
        except ValueError:
            continue

    meta = {}
    for line in lines[split + 1 : plot]:
        if not line.strip():
            continue
        parts = [x.strip() for x in line.split(",")]
        key, values = parts[0].rstrip(":"), parts[1:]
        if key.startswith("RDP Plot for event"):
            meta["title"] = key
            m = re.search(r"event\s*#?\s*(\d+)", key, re.I)
            meta["event"] = int(m.group(1)) if m else None
        elif "CI" in key:
            meta[key] = parse_pair(values)
        else:
            meta[key] = values[0] if values else None

    color_names = [x.strip() for x in lines[plot + 1].split(",")[1:] if x.strip()]
    roles = [x.strip() for x in lines[plot + 2].split(",")[1:] if x.strip()]
    header = [x.strip() for x in lines[plot + 3].split(",")]
    rows = []
    for line in lines[plot + 4 :]:
        if not line.strip():
            continue
        parts = [x.strip() for x in line.split(",")]
        if len(parts) != len(header):
            continue
        try:
            rows.append([float(x) for x in parts])
        except ValueError:
            continue
    if not rows:
        raise ValueError("No RDP plot data rows were found.")

    df = pd.DataFrame(rows, columns=header)
    series = header[1:]
    scale = max(float(df[series].to_numpy().max()), 1.0)
    colors = [COLORS.get(name.lower(), FALLBACK[i % len(FALLBACK)]) for i, name in enumerate(color_names)]

    def as_int(key, fallback=0):
        try:
            return int(float(meta.get(key) or fallback))
        except (TypeError, ValueError):
            return fallback

    meta["Maximum X-axis value"] = as_int("Maximum X-axis value", int(df[header[0]].max()))
    meta["Beginning breakpoint site"] = as_int("Beginning breakpoint site")
    meta["Ending breakpoint site"] = as_int("Ending breakpoint site")
    return {"file": path.name, "genes": genes, "meta": meta, "roles": roles, "colors": colors, "df": df, "series": series, "scale": scale}


def make_plot(data: dict) -> go.Figure:
    df, meta, series = data["df"], data["meta"], data["series"]
    xcol = df.columns[0]
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, row_heights=[0.2, 0.8], vertical_spacing=0.07,
                        subplot_titles=("Gene map", meta.get("title", "RDP event")))

    for gene in data["genes"]:
        y = FRAME_Y.get(gene["frame"], 1.5)
        fig.add_trace(go.Bar(x=[gene["end"] - gene["start"] + 1], y=[y], base=[gene["start"]], orientation="h",
                             width=0.58, marker_color="#111827" if gene["orientation"] == 1 else "#6b7280",
                             showlegend=False,
                             hovertemplate=f"Gene {gene['start']:,}–{gene['end']:,}<br>Frame {gene['frame']}<extra></extra>"), row=1, col=1)

    for i, col in enumerate(series):
        role = data["roles"][i] if i < len(data["roles"]) else "Pairwise comparison"
        fig.add_trace(go.Scatter(x=df[xcol], y=df[col] / data["scale"], mode="lines", name=role,
                                 line=dict(color=data["colors"][i % len(data["colors"])], width=2.6),
                                 customdata=df[col],
                                 hovertemplate=f"Position %{{x:,.0f}}<br>Pairwise identity %{{y:.3f}}<br>{col}<extra></extra>"), row=2, col=1)

    bands = [
        ("Beginning breakpoint 99% CI", "#818cf8", 0.12),
        ("Beginning breakpoint 95% CI", "#4f46e5", 0.18),
        ("Ending breakpoint 99% CI", "#fb7185", 0.12),
        ("Ending breakpoint 95% CI", "#e11d48", 0.18),
    ]
    for key, color, opacity in bands:
        bounds = meta.get(key)
        if isinstance(bounds, (list, tuple)) and len(bounds) == 2 and None not in bounds:
            fig.add_vrect(x0=bounds[0], x1=bounds[1], fillcolor=color, opacity=opacity, line_width=0, row=2, col=1)

    for key, color in [("Beginning breakpoint site", "#4338ca"), ("Ending breakpoint site", "#be123c")]:
        if meta.get(key):
            fig.add_vline(x=meta[key], line_color=color, line_width=2, row=2, col=1)

    fig.update_layout(template="plotly_white", height=820, hovermode="x unified", barmode="overlay",
                      margin=dict(l=60, r=30, t=80, b=55), legend=dict(orientation="h", y=1.02, x=1, xanchor="right"))
    fig.update_xaxes(range=[0, meta["Maximum X-axis value"]], title_text=meta.get("X-axis label", "Position in alignment"), row=2, col=1)
    fig.update_yaxes(tickmode="array", tickvals=[0.5, 1.5, 2.5], ticktext=["Bottom", "Middle", "Top"], range=[-0.1, 3.1], row=1, col=1)
    fig.update_yaxes(range=[0, 1.05], title_text=meta.get("Y-axis label", "Pairwise identity"), row=2, col=1)
    return fig


app = Dash(__name__, external_stylesheets=[dbc.themes.FLATLY], title="RDP Event Viewer")
server = app.server
app.layout = dbc.Container([
    dcc.Store(id="folder", data=str(DEFAULT_DIR.resolve())),
    dcc.Interval(id="poll", interval=POLL_MS),
    html.H1("RDP Event Viewer", className="display-5 fw-bold mt-4"),
    html.P("Automatically turn RDP CSV Code 1 exports into an interactive event plot.", className="lead text-secondary"),
    dbc.Card(dbc.CardBody([
        dbc.Row([
            dbc.Col([dbc.Label("RDP CSV folder"), dbc.Input(id="folder-input", value=str(DEFAULT_DIR.resolve()), debounce=True)], md=9),
            dbc.Col(dbc.Button("Use folder", id="use-folder", color="primary", className="w-100 mt-4"), md=3),
        ]),
        html.Div(id="folder-status", className="mt-2"), html.Hr(),
        dbc.Row([
            dbc.Col([dbc.Label("Event CSV"), dcc.Dropdown(id="file-picker", clearable=False)], md=8),
            dbc.Col(dbc.Checklist(id="follow-newest", options=[{"label": " Automatically open newest CSV", "value": "yes"}], value=["yes"], switch=True, className="mt-4"), md=4),
        ]),
    ]), className="mb-3"),
    html.Div(id="load-status"), html.Div(id="summary", className="my-3"),
    dcc.Graph(id="event-plot", config={"displaylogo": False, "toImageButtonOptions": {"format": "png", "filename": "rdp_event_plot", "scale": 2}}),
], fluid="xl")


@callback(Output("folder", "data"), Output("folder-status", "children"), Input("use-folder", "n_clicks"), State("folder-input", "value"), prevent_initial_call=True)
def set_folder(_, value):
    try:
        p = folder_path(value)
        return str(p), dbc.Alert(f"Watching {p}", color="success", className="py-2")
    except Exception as exc:
        return no_update, dbc.Alert(str(exc), color="danger", className="py-2")


@callback(Output("file-picker", "options"), Output("file-picker", "value"), Input("poll", "n_intervals"), Input("folder", "data"), State("file-picker", "value"), State("follow-newest", "value"))
def refresh(_, folder, selected, follow):
    try:
        files = discover(folder)
    except Exception:
        return [], None
    options = [{"label": p.name, "value": p.name} for p in files]
    names = {p.name for p in files}
    if files and ("yes" in (follow or []) or selected not in names):
        selected = files[0].name
    return options, selected


@callback(Output("event-plot", "figure"), Output("summary", "children"), Output("load-status", "children"), Input("file-picker", "value"), State("folder", "data"))
def render(filename, folder):
    if not filename:
        return go.Figure(), "", ""
    try:
        root = folder_path(folder)
        path = (root / filename).resolve()
        if path.parent != root or path.suffix.lower() != ".csv":
            raise ValueError("Invalid CSV path.")
        data = parse_code1(path)
        m = data["meta"]
        cards = dbc.Row([
            dbc.Col(dbc.Alert(f"Event #{m.get('event', '?')}", color="primary")),
            dbc.Col(dbc.Alert(f"Breakpoints {m.get('Beginning breakpoint site', '?')}–{m.get('Ending breakpoint site', '?')}", color="warning")),
            dbc.Col(dbc.Alert(f"{len(data['genes'])} genes", color="secondary")),
            dbc.Col(dbc.Alert(f"{len(data['df']):,} plotted positions", color="info")),
        ], className="g-2")
        return make_plot(data), cards, dbc.Alert(f"Loaded {filename}", color="success", className="py-2")
    except Exception as exc:
        return go.Figure(), "", dbc.Alert(str(exc), color="danger", className="py-2")


if __name__ == "__main__":
    app.run(debug=False, host="127.0.0.1", port=8050)
