from __future__ import annotations

import os
from datetime import datetime
from pathlib import Path

import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash import Dash, Input, Output, State, callback, dash_table, dcc, html, no_update

from app import empty_figure, inspect_rdp5, rdp5_alignment_overview


POLL_MS = int(os.environ.get("RDP5_POLL_MS", "3000"))
DEFAULT_PROJECT_DIR = Path(
    os.environ.get("RDP5_DIR", Path(__file__).with_name("rdp5_projects"))
).expanduser()

# The default drop-folder is created automatically. User-entered folders are
# never created implicitly; a typo should fail loudly instead of making a new
# empty directory.
if "RDP5_DIR" not in os.environ:
    DEFAULT_PROJECT_DIR.mkdir(parents=True, exist_ok=True)


def normalize_dir(value: str | os.PathLike[str]) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.exists():
        raise ValueError(f"Folder does not exist: {path}")
    if not path.is_dir():
        raise ValueError(f"Not a folder: {path}")
    return path


def scan_rdp5_folder(folder: str | os.PathLike[str]) -> list[dict]:
    root = normalize_dir(folder)
    rows = []
    for path in root.glob("*.rdp5"):
        try:
            stat = path.stat()
        except OSError:
            continue
        rows.append(
            {
                "name": path.name,
                "path": str(path.resolve()),
                "mtime_ns": stat.st_mtime_ns,
                "modified": datetime.fromtimestamp(stat.st_mtime).strftime("%Y-%m-%d %H:%M:%S"),
                "size_mb": round(stat.st_size / (1024**2), 1),
            }
        )
    rows.sort(key=lambda row: (row["mtime_ns"], row["name"].lower()), reverse=True)
    return rows


def safe_project_path(folder: str | os.PathLike[str], filename: str) -> Path:
    root = normalize_dir(folder)
    candidate = (root / filename).resolve()
    if candidate.parent != root:
        raise ValueError("Selected project is outside the configured RDP5 folder.")
    if candidate.suffix.lower() != ".rdp5":
        raise ValueError("Only .rdp5 project files are allowed.")
    if not candidate.is_file():
        raise ValueError(f"Project no longer exists: {candidate.name}")
    return candidate


def project_figures(data: dict) -> tuple[go.Figure, go.Figure, go.Figure, go.Figure]:
    overview = data.get("overview", {})
    labels = pd.DataFrame(data.get("labels", []))
    x = overview.get("x", [])

    heatmap = go.Figure(
        go.Heatmap(
            z=overview.get("matrix", []),
            x=x,
            y=overview.get("names", []),
            colorscale=[
                [0, "#f7fbff"],
                [0.15, "#c6dbef"],
                [0.35, "#6baed6"],
                [0.65, "#fdae61"],
                [1, "#b2182b"],
            ],
            zmin=0,
            zmax=35,
            colorbar=dict(title="Divergence<br>from consensus (%)"),
            hovertemplate="%{y}<br>Site %{x:,}<br>Divergence %{z:.1f}%<extra></extra>",
        )
    )
    heatmap.update_layout(
        title="Alignment divergence overview",
        xaxis_title="Alignment position (nt)",
        yaxis_title="Representative sequences",
        height=760,
    )

    profile = go.Figure()
    profile.add_trace(
        go.Scatter(
            x=x,
            y=overview.get("diversity", []),
            name="Mean divergence",
            line=dict(width=3),
            fill="tozeroy",
        )
    )
    profile.add_trace(
        go.Scatter(
            x=x,
            y=overview.get("gaps", []),
            name="Gap / uncalled sites",
            line=dict(width=2),
        )
    )
    profile.update_layout(
        title="Genome-wide signal profile",
        xaxis_title="Alignment position (nt)",
        yaxis_title="Percent",
        hovermode="x unified",
        height=500,
    )

    if labels.empty:
        subtype = empty_figure("No sequence-label metadata could be extracted")
        country = empty_figure("No country metadata could be extracted")
    else:
        subtype_counts = (
            labels.loc[labels["Subtype"].ne("Unknown"), "Subtype"]
            .value_counts()
            .head(24)
            .rename_axis("Subtype")
            .reset_index(name="Sequences")
        )
        subtype = (
            px.bar(
                subtype_counts,
                x="Subtype",
                y="Sequences",
                title="Sequence composition by subtype",
            )
            if not subtype_counts.empty
            else empty_figure("No subtype labels were detected")
        )

        country_counts = (
            labels.loc[labels["Country code"].ne("Unknown"), "Country code"]
            .value_counts()
            .head(30)
            .rename_axis("Country")
            .reset_index(name="Sequences")
        )
        country = (
            px.bar(
                country_counts,
                x="Country",
                y="Sequences",
                title="Sequence composition by country code",
            )
            if not country_counts.empty
            else empty_figure("No country codes were detected")
        )

    for figure in (heatmap, profile, subtype, country):
        figure.update_layout(
            template="plotly_white",
            margin=dict(l=50, r=25, t=70, b=50),
            hoverlabel=dict(namelength=-1),
        )
    return heatmap, profile, subtype, country


app = Dash(__name__, external_stylesheets=[dbc.themes.FLATLY], title="RDP5 Folder Visualizer")
server = app.server

app.layout = dbc.Container(
    [
        dcc.Store(id="active-folder", data=str(DEFAULT_PROJECT_DIR.resolve())),
        dcc.Store(id="project-data"),
        dcc.Interval(id="folder-poll", interval=POLL_MS, n_intervals=0),
        html.Div(
            [
                html.H1("RDP5 Folder Visualizer", className="display-5 fw-bold mb-1"),
                html.P(
                    "Point the app at an RDP output folder. New .rdp5 projects appear automatically.",
                    className="lead text-secondary",
                ),
            ],
            className="mt-4 mb-3",
        ),
        dbc.Card(
            dbc.CardBody(
                [
                    dbc.Row(
                        [
                            dbc.Col(
                                [
                                    dbc.Label("RDP5 project folder", className="fw-semibold"),
                                    dbc.Input(
                                        id="folder-input",
                                        type="text",
                                        value=str(DEFAULT_PROJECT_DIR.resolve()),
                                        debounce=True,
                                    ),
                                ],
                                md=9,
                            ),
                            dbc.Col(
                                dbc.Button(
                                    "Use folder",
                                    id="use-folder",
                                    color="primary",
                                    className="w-100 mt-4",
                                ),
                                md=3,
                            ),
                        ],
                        className="g-2",
                    ),
                    html.Div(id="folder-message", className="mt-2"),
                    html.Hr(),
                    dbc.Row(
                        [
                            dbc.Col(
                                [
                                    dbc.Label("Project", className="fw-semibold"),
                                    dcc.Dropdown(
                                        id="project-picker",
                                        placeholder="Waiting for an .rdp5 project…",
                                        clearable=False,
                                    ),
                                ],
                                md=8,
                            ),
                            dbc.Col(
                                [
                                    dbc.Label("Folder behavior", className="fw-semibold"),
                                    dbc.Checklist(
                                        id="follow-newest",
                                        options=[
                                            {
                                                "label": " Automatically open newest project",
                                                "value": "follow",
                                            }
                                        ],
                                        value=["follow"],
                                        switch=True,
                                    ),
                                ],
                                md=4,
                            ),
                        ],
                        className="g-3",
                    ),
                    dbc.Row(
                        [
                            dbc.Col(html.Div(id="scan-status", className="small text-secondary mt-2")),
                            dbc.Col(
                                dbc.Button(
                                    "Re-read selected file",
                                    id="reload-project",
                                    color="secondary",
                                    outline=True,
                                    size="sm",
                                    className="float-end mt-2",
                                )
                            ),
                        ]
                    ),
                ]
            ),
            className="mb-3",
        ),
        html.Div(id="project-status"),
        html.Div(id="summary-cards", className="my-3"),
        dbc.Tabs(
            [
                dbc.Tab(
                    [
                        dcc.Graph(
                            id="alignment-heatmap",
                            config={
                                "displaylogo": False,
                                "toImageButtonOptions": {
                                    "format": "png",
                                    "filename": "rdp5_alignment_overview",
                                    "scale": 2,
                                },
                            },
                        ),
                        dcc.Graph(
                            id="genome-profile",
                            config={
                                "displaylogo": False,
                                "toImageButtonOptions": {
                                    "format": "png",
                                    "filename": "rdp5_genome_profile",
                                    "scale": 2,
                                },
                            },
                        ),
                    ],
                    label="Alignment overview",
                ),
                dbc.Tab(
                    [
                        dbc.Row(
                            [
                                dbc.Col(dcc.Graph(id="subtype-chart"), lg=6),
                                dbc.Col(dcc.Graph(id="country-chart"), lg=6),
                            ],
                            className="g-2 mt-1",
                        )
                    ],
                    label="Dataset composition",
                ),
                dbc.Tab(
                    [
                        dash_table.DataTable(
                            id="label-table",
                            page_size=20,
                            sort_action="native",
                            filter_action="native",
                            export_format="csv",
                            style_table={"overflowX": "auto"},
                            style_cell={
                                "fontFamily": "system-ui",
                                "fontSize": 13,
                                "padding": "7px",
                                "textAlign": "left",
                            },
                        )
                    ],
                    label="Sequence labels",
                ),
                dbc.Tab(
                    dbc.Card(
                        dbc.CardBody(
                            [
                                html.H4("How this mode works"),
                                dcc.Markdown(
                                    """
1. RDP writes or saves an `.rdp5` project into the configured directory.
2. This page rescans that directory every few seconds.
3. With **Automatically open newest project** enabled, the newest project is selected as soon as it appears.
4. The server reads the binary RDP5 file directly from disk and extracts the embedded sequence labels and alignment for visualization.

The browser never needs to upload the RDP5 file.

**Current limitation:** the RDP5 binary format is not publicly documented. This app reads the structures we can validate reliably from the example project—sequence labels and the embedded alignment. It does not yet claim to decode RDP's saved event calls from the binary project.
"""
                                ),
                            ]
                        ),
                        className="my-3",
                    ),
                    label="About",
                ),
            ]
        ),
        html.P(
            "Local teaching/exploration tool. Bind remains on 127.0.0.1 so the folder browser is not exposed to the network.",
            className="small text-secondary text-center my-4",
        ),
    ],
    fluid="xl",
)


@callback(
    Output("active-folder", "data"),
    Output("folder-message", "children"),
    Input("use-folder", "n_clicks"),
    State("folder-input", "value"),
    prevent_initial_call=True,
)
def choose_folder(_, folder_value):
    try:
        folder = normalize_dir(folder_value or "")
        return str(folder), dbc.Alert(
            f"Watching {folder}", color="success", className="py-2 mb-0"
        )
    except Exception as exc:
        return no_update, dbc.Alert(str(exc), color="danger", className="py-2 mb-0")


@callback(
    Output("project-picker", "options"),
    Output("project-picker", "value"),
    Output("scan-status", "children"),
    Input("folder-poll", "n_intervals"),
    Input("active-folder", "data"),
    State("project-picker", "value"),
    State("follow-newest", "value"),
)
def discover_projects(_, folder, selected, follow):
    try:
        projects = scan_rdp5_folder(folder)
    except Exception as exc:
        return [], None, f"Folder scan failed: {exc}"
    if not projects:
        return [], None, f"No .rdp5 files found in {folder}"

    options = [
        {
            "label": f'{row["name"]}  ·  {row["size_mb"]:.1f} MB  ·  {row["modified"]}',
            "value": row["name"],
        }
        for row in projects
    ]
    names = {row["name"] for row in projects}
    if "follow" in (follow or []) or selected not in names:
        selected = projects[0]["name"]

    return options, selected, f'{len(projects)} project(s) found · watching every {POLL_MS / 1000:g} s'


@callback(
    Output("project-data", "data"),
    Output("project-status", "children"),
    Input("project-picker", "value"),
    Input("reload-project", "n_clicks"),
    State("active-folder", "data"),
    prevent_initial_call=True,
)
def load_project(filename, _, folder):
    if not filename:
        return no_update, ""
    try:
        path = safe_project_path(folder, filename)
        raw = path.read_bytes()
        labels = inspect_rdp5(raw)
        overview = rdp5_alignment_overview(raw, labels)
        # Keep the browser payload small. The folder view does not need to ship
        # all embedded sequences to the client.
        overview.pop("sequences_raw", None)
        stat = path.stat()
        data = {
            "filename": path.name,
            "path": str(path),
            "bytes": len(raw),
            "modified": datetime.fromtimestamp(stat.st_mtime).isoformat(timespec="seconds"),
            "labels": labels.to_dict("records"),
            "overview": overview,
        }
        message = (
            f'Loaded {path.name}: {overview.get("sequences", 0):,} sequences × '
            f'{overview.get("length", 0):,} sites.'
        )
        return data, dbc.Alert(message, color="success", className="py-2")
    except Exception as exc:
        return no_update, dbc.Alert(str(exc), color="danger", className="py-2")


@callback(
    Output("alignment-heatmap", "figure"),
    Output("genome-profile", "figure"),
    Output("subtype-chart", "figure"),
    Output("country-chart", "figure"),
    Output("label-table", "data"),
    Output("label-table", "columns"),
    Output("summary-cards", "children"),
    Input("project-data", "data"),
)
def render_project(data):
    blank = empty_figure("Waiting for an RDP5 project")
    if not data:
        return blank, blank, blank, blank, [], [], ""

    heatmap, profile, subtype, country = project_figures(data)
    frame = pd.DataFrame(data.get("labels", []))
    overview = data.get("overview", {})
    known_subtypes = (
        frame["Subtype"].replace("Unknown", np.nan).nunique() if not frame.empty else 0
    )
    known_countries = (
        frame["Country code"].replace("Unknown", np.nan).nunique() if not frame.empty else 0
    )
    summary = dbc.Row(
        [
            dbc.Col(dbc.Alert(f'{overview.get("sequences", 0):,} sequences', color="primary")),
            dbc.Col(dbc.Alert(f'{overview.get("length", 0):,} sites', color="secondary")),
            dbc.Col(dbc.Alert(f"{known_subtypes:,} subtypes", color="info")),
            dbc.Col(dbc.Alert(f"{known_countries:,} country codes", color="success")),
        ],
        className="g-2",
    )
    return (
        heatmap,
        profile,
        subtype,
        country,
        frame.to_dict("records"),
        [{"name": column, "id": column} for column in frame.columns],
        summary,
    )


if __name__ == "__main__":
    app.run(debug=True, host="127.0.0.1", port=8050)
