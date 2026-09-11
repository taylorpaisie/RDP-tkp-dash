from __future__ import annotations

import base64
import json
import math
import os
from typing import Any

import dash_bootstrap_components as dbc
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash import Dash, Input, Output, State, callback, dash_table, dcc, html, no_update


def decode_json(contents: str) -> dict[str, Any]:
    raw = base64.b64decode(contents.split(",", 1)[1])
    result = json.loads(raw.decode("utf-8-sig"))
    if not isinstance(result, dict):
        raise ValueError("The Datamonkey file must contain a JSON object at its top level.")
    return result


def method_name(doc: dict[str, Any]) -> str:
    analysis = doc.get("analysis", {})
    candidates = [analysis.get("info"), analysis.get("name"), doc.get("method"), doc.get("analysis")]
    for value in candidates:
        if isinstance(value, str) and value.strip():
            return value.strip()
        if isinstance(value, dict):
            for key in ("method", "name"):
                if isinstance(value.get(key), str):
                    return value[key]
    text = json.dumps(doc.get("test results", {})).lower()
    for known in ("busted", "meme", "fel", "fubar", "slac", "relax", "absrel"):
        if known in text:
            return known.upper()
    return "HyPhy / Datamonkey"


def input_summary(doc: dict[str, Any]) -> dict[str, Any]:
    inp = doc.get("input", {}) if isinstance(doc.get("input"), dict) else {}
    return {
        "sequences": inp.get("number of sequences") or inp.get("sequences") or "—",
        "sites": inp.get("number of sites") or inp.get("sites") or "—",
        "partitions": inp.get("partition count") or inp.get("partitions") or len(doc.get("data partitions", {}) or {}) or "—",
    }


def _headers_and_rows(section: Any) -> tuple[list[str], list[list[Any]]]:
    if not isinstance(section, dict):
        return [], []
    headers = section.get("headers", [])
    content = section.get("content", {})
    names = []
    for h in headers if isinstance(headers, list) else []:
        names.append(str(h[0] if isinstance(h, list) and h else h))
    rows: list[list[Any]] = []
    if isinstance(content, dict):
        for partition, values in content.items():
            if isinstance(values, list):
                for i, row in enumerate(values, 1):
                    if isinstance(row, list):
                        rows.append([str(partition), i, *row])
    elif isinstance(content, list):
        for i, row in enumerate(content, 1):
            if isinstance(row, list):
                rows.append(["1", i, *row])
    return ["Partition", "Site", *names], rows


def site_table(doc: dict[str, Any]) -> pd.DataFrame:
    for key in ("MLE", "site results", "site-results"):
        section = doc.get(key)
        if isinstance(section, dict) and "headers" in section:
            cols, rows = _headers_and_rows(section)
            if rows:
                width = min(len(cols), min(map(len, rows)))
                return pd.DataFrame([r[:width] for r in rows], columns=cols[:width])
    mle = doc.get("MLE", {})
    if isinstance(mle, dict):
        for value in mle.values():
            cols, rows = _headers_and_rows(value)
            if rows:
                width = min(len(cols), min(map(len, rows)))
                return pd.DataFrame([r[:width] for r in rows], columns=cols[:width])
    return pd.DataFrame()


def model_table(doc: dict[str, Any]) -> pd.DataFrame:
    fits = doc.get("fits", {})
    rows = []
    if isinstance(fits, dict):
        for model, values in fits.items():
            row = {"Model": model}
            if isinstance(values, dict):
                for key, value in values.items():
                    if isinstance(value, (str, int, float, bool)) or value is None:
                        row[str(key)] = value
            rows.append(row)
    return pd.DataFrame(rows)


def branch_table(doc: dict[str, Any]) -> pd.DataFrame:
    attrs = doc.get("branch attributes", {})
    rows = []
    if isinstance(attrs, dict):
        for partition, branches in attrs.items():
            if not isinstance(branches, dict):
                continue
            for branch, values in branches.items():
                row = {"Partition": partition, "Branch": branch}
                if isinstance(values, dict):
                    for key, value in values.items():
                        if isinstance(value, (str, int, float, bool)) or value is None:
                            row[str(key)] = value
                rows.append(row)
    return pd.DataFrame(rows)


def test_rows(doc: dict[str, Any]) -> list[dict[str, Any]]:
    tests = doc.get("test results", {})
    if not isinstance(tests, dict):
        return []
    rows = []
    for key, value in tests.items():
        if isinstance(value, (str, int, float, bool)) or value is None:
            rows.append({"Statistic": key, "Value": value})
    return rows


def find_column(frame: pd.DataFrame, terms: tuple[str, ...]) -> str | None:
    for col in frame.columns:
        normalized = str(col).lower().replace("_", " ").replace("-", " ")
        if any(term in normalized for term in terms):
            return str(col)
    return None


def site_figure(frame: pd.DataFrame) -> go.Figure:
    if frame.empty:
        return go.Figure().update_layout(template="plotly_white", title="No site-level table found in this result")
    p_col = find_column(frame, ("p value", "p-value", "q value", "posterior probability"))
    omega_col = find_column(frame, ("omega", "dnds", "dN/dS".lower(), "beta-alpha", "beta minus alpha"))
    color_col = p_col or omega_col
    plot = frame.copy()
    if color_col:
        plot[color_col] = pd.to_numeric(plot[color_col], errors="coerce")
    if p_col:
        plot["Evidence"] = -plot[p_col].clip(lower=1e-300).map(math.log10)
        fig = px.scatter(plot, x="Site", y="Evidence", color="Partition", hover_data=list(plot.columns),
                         title=f"Site-level evidence (−log10 {p_col})")
        fig.add_hline(y=-math.log10(0.05), line_dash="dash", line_color="#dc2626", annotation_text="0.05")
    elif omega_col:
        fig = px.scatter(plot, x="Site", y=omega_col, color="Partition", hover_data=list(plot.columns),
                         title=f"Site-level selection signal ({omega_col})")
        fig.add_hline(y=1, line_dash="dash", line_color="#dc2626", annotation_text="neutral")
    else:
        numeric = [c for c in plot.columns if c not in ("Partition", "Site") and pd.to_numeric(plot[c], errors="coerce").notna().any()]
        if not numeric:
            return go.Figure().update_layout(template="plotly_white", title="Site table loaded; choose the table tab to inspect it")
        plot[numeric[0]] = pd.to_numeric(plot[numeric[0]], errors="coerce")
        fig = px.scatter(plot, x="Site", y=numeric[0], color="Partition", hover_data=list(plot.columns), title=f"Site-level results: {numeric[0]}")
    fig.update_layout(template="plotly_white", height=520, margin=dict(l=55, r=25, t=70, b=50))
    return fig


app = Dash(__name__, external_stylesheets=[dbc.themes.FLATLY], title="Datamonkey JSON Visualizer")
server = app.server
app.layout = dbc.Container([
    dcc.Store(id="dm-data"),
    html.H1("Datamonkey JSON Visualizer", className="display-5 fw-bold mt-4"),
    html.P("Turn HyPhy/Datamonkey result JSON into readable evidence plots and tables.", className="lead text-secondary"),
    dbc.Card(dbc.CardBody([
        dcc.Upload(id="dm-upload", accept=".json,application/json", children=html.Div(["Drop a Datamonkey JSON file here — or ", html.A("browse")]), className="border rounded p-5 text-center"),
        html.Div(id="dm-status", className="mt-3"),
    ]), className="mb-3"),
    html.Div(id="dm-cards", className="mb-3"),
    dbc.Tabs([
        dbc.Tab(dcc.Graph(id="dm-sites", config={"displaylogo": False}), label="Selection by site"),
        dbc.Tab([html.H5("Test results", className="mt-3"), dash_table.DataTable(id="dm-tests", style_table={"overflowX": "auto"}, style_cell={"padding": "8px", "textAlign": "left"}),
                 html.H5("Model fits", className="mt-4"), dash_table.DataTable(id="dm-fits", page_size=12, sort_action="native", style_table={"overflowX": "auto"}, style_cell={"padding": "7px", "fontSize": 13})], label="Global tests & models"),
        dbc.Tab(dash_table.DataTable(id="dm-site-table", page_size=18, sort_action="native", filter_action="native", export_format="csv", style_table={"overflowX": "auto", "marginTop": "1rem"}, style_cell={"padding": "7px", "fontSize": 13}), label="Site table"),
        dbc.Tab(dash_table.DataTable(id="dm-branches", page_size=18, sort_action="native", filter_action="native", export_format="csv", style_table={"overflowX": "auto", "marginTop": "1rem"}, style_cell={"padding": "7px", "fontSize": 13}), label="Branches"),
        dbc.Tab(html.Pre(id="dm-raw", style={"maxHeight": "650px", "overflow": "auto", "fontSize": "12px", "marginTop": "1rem"}), label="JSON explorer"),
    ]),
    dbc.Alert("Interpret thresholds according to the selected Datamonkey method; a highlighted site is evidence, not a tiny evolutionary confession.", color="warning", className="mt-3"),
], fluid="xl")


@callback(Output("dm-data", "data"), Output("dm-status", "children"), Input("dm-upload", "contents"), State("dm-upload", "filename"), prevent_initial_call=True)
def ingest(contents, filename):
    if not contents:
        return no_update, no_update
    try:
        doc = decode_json(contents)
        payload = {"filename": filename, "method": method_name(doc), "document": doc}
        return payload, dbc.Alert(f"Loaded {filename or 'JSON'} · detected {payload['method']}", color="success", className="py-2")
    except Exception as exc:
        return no_update, dbc.Alert(f"Could not read this JSON: {exc}", color="danger", className="py-2")


@callback(Output("dm-cards", "children"), Output("dm-sites", "figure"), Output("dm-tests", "data"), Output("dm-tests", "columns"),
          Output("dm-fits", "data"), Output("dm-fits", "columns"), Output("dm-site-table", "data"), Output("dm-site-table", "columns"),
          Output("dm-branches", "data"), Output("dm-branches", "columns"), Output("dm-raw", "children"), Input("dm-data", "data"))
def render(data):
    blank = go.Figure().update_layout(template="plotly_white", title="Upload a Datamonkey JSON result")
    if not data:
        return "", blank, [], [], [], [], [], [], [], [], ""
    doc = data["document"]
    sites, fits, branches = site_table(doc), model_table(doc), branch_table(doc)
    tests = test_rows(doc)
    summary = input_summary(doc)
    cards = dbc.Row([
        dbc.Col(dbc.Alert(data["method"], color="primary")),
        dbc.Col(dbc.Alert(f"{summary['sequences']} sequences", color="secondary")),
        dbc.Col(dbc.Alert(f"{summary['sites']} sites", color="info")),
        dbc.Col(dbc.Alert(f"{len(sites):,} site rows", color="success")),
    ])
    def cols(frame): return [{"name": str(c), "id": str(c)} for c in frame.columns]
    return (cards, site_figure(sites), tests, [{"name": "Statistic", "id": "Statistic"}, {"name": "Value", "id": "Value"}],
            fits.to_dict("records"), cols(fits), sites.to_dict("records"), cols(sites), branches.to_dict("records"), cols(branches),
            json.dumps(doc, indent=2, ensure_ascii=False)[:2_000_000])


if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=int(os.environ.get("PORT", "8051")))
