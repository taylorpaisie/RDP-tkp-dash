from __future__ import annotations

import base64
import io
import re
from pathlib import Path

import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash import Dash, Input, Output, State, callback, dash_table, dcc, html, no_update


METHOD_NAMES = ["RDP", "GENECONV", "BootScan", "MaxChi", "Chimaera", "SiScan", "3Seq", "LARD", "PhylPro"]
COLUMN_ALIASES = {
    "event": ["event", "eventnumber", "eventno", "event#", "recombinationevent"],
    "recombinant": ["recombinant", "recombinantsequence", "recseq", "sequence"],
    "major_parent": ["majorparent", "major_parent", "parent1", "parent_1", "major"],
    "minor_parent": ["minorparent", "minor_parent", "parent2", "parent_2", "minor"],
    "start": ["beginningbreakpoint", "beginbreakpoint", "startbreakpoint", "breakpointstart", "start", "begin"],
    "end": ["endingbreakpoint", "endbreakpoint", "stopbreakpoint", "breakpointend", "end", "stop"],
}


def _key(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def _decode_upload(contents: str) -> bytes:
    return base64.b64decode(contents.split(",", 1)[1])


def _empty_figure(title: str) -> go.Figure:
    fig = go.Figure()
    fig.update_layout(
        template="plotly_white",
        title=title,
        height=420,
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        annotations=[dict(text=title, x=.5, y=.5, xref="paper", yref="paper", showarrow=False,
                          font=dict(size=17, color="#64748b"))],
        margin=dict(l=35, r=35, t=65, b=35),
    )
    return fig


def _read_table(raw: bytes) -> pd.DataFrame:
    text = raw.decode("utf-8-sig", errors="replace")
    attempts = [
        dict(sep=None, engine="python"),
        dict(sep=",", engine="python"),
        dict(sep="\t", engine="python"),
        dict(sep=";", engine="python"),
    ]
    last = None
    for kwargs in attempts:
        try:
            frame = pd.read_csv(io.StringIO(text), **kwargs)
            if frame.shape[1] >= 3:
                return frame
        except Exception as exc:
            last = exc
    raise ValueError(f"Could not read the file as an RDP result table: {last}")


def _find_column(columns: list[str], aliases: list[str]) -> str | None:
    normalized = {_key(c): c for c in columns}
    for alias in aliases:
        if _key(alias) in normalized:
            return normalized[_key(alias)]
    return None


def _method_support(value) -> bool:
    if pd.isna(value):
        return False
    text = str(value).strip().lower()
    if text in {"", "-", "--", "na", "n/a", "none", "false", "0"}:
        return False
    number = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.notna(number):
        # RDP method fields commonly contain p-values; any finite positive value is evidence that method reported the event.
        return bool(number > 0)
    return True


def parse_rdp_csv(raw: bytes) -> tuple[pd.DataFrame, list[str]]:
    frame = _read_table(raw)
    original_columns = [str(c) for c in frame.columns]

    rename = {}
    for target, aliases in COLUMN_ALIASES.items():
        found = _find_column(original_columns, aliases)
        if found:
            rename[found] = target
    frame = frame.rename(columns=rename)

    required = {"recombinant", "start", "end"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(
            "This does not look like an RDP recombination-event CSV. "
            f"I could not identify: {', '.join(sorted(missing))}. "
            "Expected columns equivalent to Recombinant, Beginning breakpoint, and Ending breakpoint."
        )

    for column in ["start", "end"]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame = frame.dropna(subset=["start", "end", "recombinant"]).copy()
    if frame.empty:
        raise ValueError("No rows contained usable recombinant and breakpoint coordinates.")

    frame["start"] = frame["start"].astype(int)
    frame["end"] = frame["end"].astype(int)
    frame["left"] = frame[["start", "end"]].min(axis=1)
    frame["right"] = frame[["start", "end"]].max(axis=1)
    frame["length"] = frame["right"] - frame["left"] + 1
    frame["recombinant"] = frame["recombinant"].astype(str).str.strip()

    if "event" not in frame:
        frame["event"] = np.arange(1, len(frame) + 1)
    if "major_parent" not in frame:
        frame["major_parent"] = "Unknown"
    if "minor_parent" not in frame:
        frame["minor_parent"] = "Unknown"
    frame["major_parent"] = frame["major_parent"].fillna("Unknown").astype(str).str.strip()
    frame["minor_parent"] = frame["minor_parent"].fillna("Unknown").astype(str).str.strip()

    methods = []
    method_lookup = {_key(c): c for c in frame.columns}
    for canonical in METHOD_NAMES:
        candidates = [canonical, f"{canonical}pvalue", f"{canonical}p-value", f"{canonical}p"]
        found = next((method_lookup.get(_key(x)) for x in candidates if _key(x) in method_lookup), None)
        if found and found not in methods:
            methods.append(found)

    if not methods:
        # Fallback: method-like columns that are mostly blank / numeric p-values / short flags.
        protected = {"event", "recombinant", "major_parent", "minor_parent", "start", "end", "left", "right", "length"}
        for column in frame.columns:
            if column in protected:
                continue
            sample = frame[column].dropna().astype(str).head(30)
            if not sample.empty and sample.str.len().median() <= 18:
                numeric_fraction = pd.to_numeric(sample, errors="coerce").notna().mean()
                if numeric_fraction >= .5:
                    methods.append(column)

    for method in methods:
        frame[f"support::{method}"] = frame[method].map(_method_support)
    support_cols = [c for c in frame.columns if c.startswith("support::")]
    frame["method_count"] = frame[support_cols].sum(axis=1) if support_cols else 0

    return frame, methods


def event_map(frame: pd.DataFrame) -> go.Figure:
    ordered = frame.sort_values(["recombinant", "left", "right"]).copy()
    ordered["track"] = ordered["recombinant"].astype(str)
    hover = ["event", "major_parent", "minor_parent", "left", "right", "length", "method_count"]
    fig = px.bar(
        ordered,
        x="length",
        y="track",
        base="left",
        color="minor_parent",
        orientation="h",
        hover_data=hover,
        title="Recombination event map",
        labels={"track": "Recombinant", "minor_parent": "Minor parent", "length": "Event span (nt)"},
    )
    fig.update_traces(marker_line_width=.5, marker_line_color="white")
    fig.update_layout(barmode="overlay", height=max(480, 34 * ordered["track"].nunique() + 150))
    fig.update_xaxes(title="Genome / alignment position (nt)")
    return fig


def breakpoint_density(frame: pd.DataFrame) -> go.Figure:
    points = pd.concat([
        frame[["left"]].rename(columns={"left": "position"}).assign(kind="Start breakpoint"),
        frame[["right"]].rename(columns={"right": "position"}).assign(kind="End breakpoint"),
    ], ignore_index=True)
    bins = min(80, max(20, int(np.sqrt(len(points)) * 3)))
    fig = px.histogram(points, x="position", color="kind", nbins=bins, barmode="overlay",
                       title="Breakpoint density", labels={"position": "Genome / alignment position (nt)", "count": "Breakpoints"})
    fig.update_traces(opacity=.72)
    return fig


def method_heatmap(frame: pd.DataFrame, methods: list[str]) -> go.Figure:
    if not methods:
        return _empty_figure("No recognizable per-method support columns were found")
    support = pd.DataFrame({method: frame[f"support::{method}"].astype(int) for method in methods})
    labels = [f"Event {row.event}: {row.recombinant}" for row in frame[["event", "recombinant"]].itertuples(index=False)]
    fig = go.Figure(go.Heatmap(
        z=support.T.values,
        x=labels,
        y=methods,
        zmin=0,
        zmax=1,
        colorscale=[[0, "#eef2f7"], [.499, "#eef2f7"], [.5, "#1f77b4"], [1, "#1f77b4"]],
        showscale=False,
        hovertemplate="%{y}<br>%{x}<br>Support: %{z}<extra></extra>",
    ))
    fig.update_layout(title="Method support matrix", height=max(400, 42 * len(methods) + 140))
    fig.update_xaxes(tickangle=-45, showticklabels=len(labels) <= 45)
    return fig


def parental_flow(frame: pd.DataFrame) -> go.Figure:
    links = (frame.groupby(["major_parent", "minor_parent", "recombinant"], dropna=False)
             .size().reset_index(name="events").sort_values("events", ascending=False))
    if links.empty:
        return _empty_figure("No parental relationships available")

    major = [f"Major · {x}" for x in links["major_parent"].astype(str).unique()]
    minor = [f"Minor · {x}" for x in links["minor_parent"].astype(str).unique()]
    rec = [f"Recombinant · {x}" for x in links["recombinant"].astype(str).unique()]
    nodes = major + minor + rec
    index = {name: i for i, name in enumerate(nodes)}

    source, target, value = [], [], []
    first = links.groupby(["major_parent", "minor_parent"], dropna=False)["events"].sum().reset_index()
    for row in first.itertuples(index=False):
        source.append(index[f"Major · {row.major_parent}"])
        target.append(index[f"Minor · {row.minor_parent}"])
        value.append(int(row.events))
    for row in links.itertuples(index=False):
        source.append(index[f"Minor · {row.minor_parent}"])
        target.append(index[f"Recombinant · {row.recombinant}"])
        value.append(int(row.events))

    fig = go.Figure(go.Sankey(
        arrangement="snap",
        node=dict(label=nodes, pad=18, thickness=18),
        link=dict(source=source, target=target, value=value),
    ))
    fig.update_layout(title="Parent → recombinant relationship flow", height=650, font=dict(size=12))
    return fig


app = Dash(__name__, external_stylesheets=[dbc.themes.FLATLY], title="RDP CSV Visualizer")
server = app.server

app.layout = dbc.Container([
    dcc.Store(id="rdp-store"),
    html.Div([
        html.H1("RDP CSV → Pretty Pictures", className="display-5 fw-bold mb-1"),
        html.P("Drop in an RDP-exported CSV/TSV and turn the event table into interactive, presentation-ready figures.",
               className="lead text-secondary"),
    ], className="mt-4 mb-3"),

    dbc.Card(dbc.CardBody([
        dcc.Upload(
            id="rdp-upload",
            children=html.Div([html.Strong("Drop an RDP CSV/TSV here"), " — or ", html.A("browse")]),
            className="border rounded p-5 text-center bg-light",
            multiple=False,
        ),
        html.Div(id="rdp-status", className="mt-3"),
    ]), className="shadow-sm mb-4"),

    html.Div(id="summary-row", className="mb-3"),

    dbc.Tabs([
        dbc.Tab([dcc.Graph(id="event-map", config={"displaylogo": False, "toImageButtonOptions": {"format": "png", "scale": 2}})],
                label="Event map"),
        dbc.Tab([dcc.Graph(id="flow-chart", config={"displaylogo": False, "toImageButtonOptions": {"format": "png", "scale": 2}})],
                label="Parent flow"),
        dbc.Tab([dcc.Graph(id="breakpoint-chart", config={"displaylogo": False, "toImageButtonOptions": {"format": "png", "scale": 2}})],
                label="Breakpoint density"),
        dbc.Tab([dcc.Graph(id="method-chart", config={"displaylogo": False, "toImageButtonOptions": {"format": "png", "scale": 2}})],
                label="Method support"),
        dbc.Tab([
            html.Div(className="my-3"),
            dash_table.DataTable(
                id="event-table",
                page_size=20,
                sort_action="native",
                filter_action="native",
                export_format="csv",
                style_table={"overflowX": "auto"},
                style_cell={"fontFamily": "system-ui", "fontSize": 12, "padding": "7px", "textAlign": "left"},
                style_header={"fontWeight": "700"},
            )
        ], label="Clean table"),
    ]),

    dbc.Alert([
        html.Strong("Tip: "),
        "Use the camera button in each Plotly toolbar to save a high-resolution PNG. This app visualizes RDP output; it does not rerun recombination detection."
    ], color="info", className="my-4"),
], fluid="xl")


@callback(
    Output("rdp-store", "data"),
    Output("rdp-status", "children"),
    Input("rdp-upload", "contents"),
    State("rdp-upload", "filename"),
    prevent_initial_call=True,
)
def ingest(contents, filename):
    if not contents:
        return no_update, no_update
    try:
        frame, methods = parse_rdp_csv(_decode_upload(contents))
        payload = {
            "filename": filename or "RDP export",
            "records": frame.to_dict("records"),
            "methods": methods,
        }
        message = (
            f"Loaded {len(frame):,} recombination events from {filename or 'RDP export'} "
            f"across {frame['recombinant'].nunique():,} recombinant sequence(s)."
        )
        return payload, dbc.Alert(message, color="success", className="py-2")
    except Exception as exc:
        return no_update, dbc.Alert(str(exc), color="danger", className="py-2")


@callback(
    Output("event-map", "figure"),
    Output("flow-chart", "figure"),
    Output("breakpoint-chart", "figure"),
    Output("method-chart", "figure"),
    Output("event-table", "data"),
    Output("event-table", "columns"),
    Output("summary-row", "children"),
    Input("rdp-store", "data"),
)
def render(data):
    blank = _empty_figure("Upload an RDP CSV/TSV to begin")
    if not data:
        return blank, blank, blank, blank, [], [], ""

    frame = pd.DataFrame(data["records"])
    methods = data.get("methods", [])
    cleaned_columns = [
        "event", "recombinant", "major_parent", "minor_parent",
        "start", "end", "length", "method_count"
    ]
    cleaned_columns = [c for c in cleaned_columns if c in frame.columns]
    cleaned = frame[cleaned_columns].copy()

    max_coord = int(frame["right"].max())
    median_span = int(frame["length"].median())
    cards = dbc.Row([
        dbc.Col(dbc.Card(dbc.CardBody([html.Div("Events", className="text-secondary small"), html.H3(f"{len(frame):,}")])), md=3),
        dbc.Col(dbc.Card(dbc.CardBody([html.Div("Recombinants", className="text-secondary small"), html.H3(f"{frame['recombinant'].nunique():,}")])), md=3),
        dbc.Col(dbc.Card(dbc.CardBody([html.Div("Genome span", className="text-secondary small"), html.H3(f"{max_coord:,} nt")])), md=3),
        dbc.Col(dbc.Card(dbc.CardBody([html.Div("Median event", className="text-secondary small"), html.H3(f"{median_span:,} nt")])), md=3),
    ], className="g-3")

    figures = [event_map(frame), parental_flow(frame), breakpoint_density(frame), method_heatmap(frame, methods)]
    for fig in figures:
        fig.update_layout(template="plotly_white", margin=dict(l=50, r=30, t=75, b=55), hoverlabel=dict(namelength=-1))

    return (
        figures[0], figures[1], figures[2], figures[3],
        cleaned.to_dict("records"), [{"name": c.replace("_", " ").title(), "id": c} for c in cleaned.columns], cards
    )


if __name__ == "__main__":
    app.run(debug=True, host="127.0.0.1", port=8050)
