from __future__ import annotations

import base64
import io
import math
import re
from collections import Counter

import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash import Dash, Input, Output, State, callback, dash_table, dcc, html, no_update


EVENT_ALIASES = {
    "method": ["method", "program"], "start": ["start", "begin", "breakpoint_start"],
    "end": ["end", "stop", "breakpoint_end"], "recombinant": ["recombinant", "recseq", "sequence"],
    "parent1": ["parent1", "major_parent", "parent_1"], "parent2": ["parent2", "minor_parent", "parent_2"],
    "pvalue": ["pvalue", "p-value", "p_value", "probability"],
}


def decode_upload(contents: str) -> bytes:
    return base64.b64decode(contents.split(",", 1)[1])


def parse_fasta(raw: bytes) -> tuple[list[str], list[str]]:
    text = raw.decode("utf-8-sig", errors="strict")
    names, seqs, current = [], [], []
    for original in text.splitlines():
        line = original.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current:
                seqs.append("".join(current).upper()); current = []
            name = line[1:].strip()
            if not name:
                raise ValueError("A FASTA record has an empty name.")
            names.append(name)
        elif not names:
            raise ValueError("Sequence data occurs before the first FASTA header.")
        else:
            current.append(re.sub(r"\s+", "", line))
    if current:
        seqs.append("".join(current).upper())
    if len(names) != len(seqs) or not seqs:
        raise ValueError("The FASTA file contains an incomplete record.")
    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        raise ValueError(f"Sequences are not aligned; observed lengths: {sorted(lengths)[:8]}")
    return names, seqs


def alignment_summary(names: list[str], seqs: list[str], window: int = 100) -> tuple[pd.DataFrame, pd.DataFrame]:
    arr = np.array([list(s) for s in seqs], dtype="U1")
    rows = []
    for name, row in zip(names, arr):
        gaps = np.isin(row, ["-", ".", "?"]).sum()
        ambiguous = (~np.isin(row, ["A", "C", "G", "T", "-", ".", "?"])).sum()
        rows.append({"Sequence": name, "Sites": len(row), "Gaps": int(gaps), "Gap %": round(100*gaps/len(row), 2),
                     "Ambiguous": int(ambiguous), "Ambiguous %": round(100*ambiguous/len(row), 2)})
    profile=[]
    for left in range(0, arr.shape[1], window):
        block=arr[:, left:left+window]
        variable=0; entropy=[]; gaps=[]
        for col in block.T:
            called=[x for x in col if x in "ACGT"]
            counts=Counter(called)
            variable += len(counts)>1
            if called:
                ps=np.array(list(counts.values()),dtype=float)/len(called)
                entropy.append(float(-(ps*np.log2(ps)).sum()))
            gaps.append(np.isin(col,["-",".","?"]).mean())
        profile.append({"Start":left+1,"End":min(left+window,arr.shape[1]),"Midpoint":left+len(block.T)/2,
                        "Variable sites":variable,"Mean entropy":np.mean(entropy) if entropy else 0,
                        "Gap %":100*np.mean(gaps) if gaps else 0})
    return pd.DataFrame(rows), pd.DataFrame(profile)


def normalize_events(raw: bytes) -> pd.DataFrame:
    text=raw.decode("utf-8-sig",errors="strict")
    frame=pd.read_csv(io.StringIO(text),sep=None,engine="python")
    normalized={re.sub(r"[^a-z0-9_-]","",str(c).lower().strip()):c for c in frame.columns}
    rename={}
    for target,aliases in EVENT_ALIASES.items():
        for alias in aliases:
            key=re.sub(r"[^a-z0-9_-]","",alias.lower())
            if key in normalized:
                rename[normalized[key]]=target.title() if target!="pvalue" else "Pvalue"; break
    frame=frame.rename(columns=rename)
    required={"Start","End","Recombinant"}
    if not required.issubset(frame.columns):
        raise ValueError("Result table needs Start, End, and Recombinant columns (common name variants are accepted).")
    frame["Start"]=pd.to_numeric(frame["Start"],errors="coerce")
    frame["End"]=pd.to_numeric(frame["End"],errors="coerce")
    if "Pvalue" in frame: frame["Pvalue"]=pd.to_numeric(frame["Pvalue"],errors="coerce")
    if "Method" not in frame: frame["Method"]="RDP event"
    return frame.dropna(subset=["Start","End"])


def inspect_rdp5(raw: bytes) -> pd.DataFrame:
    if not raw.startswith(b"RDP5 Project File"):
        raise ValueError("The file does not contain the expected RDP5 project signature.")
    chunks=re.findall(rb"[A-Za-z0-9][A-Za-z0-9_.-]{5,95}",raw)
    seen=[]
    for chunk in chunks:
        value=chunk.decode("ascii",errors="ignore").strip("._-")
        if ("." in value or "_" in value) and value not in seen and not value.isdigit(): seen.append(value)
    return pd.DataFrame({"Probable sequence/project label":seen[:2000]})


app=Dash(__name__,external_stylesheets=[dbc.themes.FLATLY],title="RDP-tkp Visualizer")
server=app.server
app.layout=dbc.Container([
    dcc.Store(id="parsed"),
    html.H1("RDP-tkp Visualizer",className="display-5 fw-bold mt-4"),
    html.P("Explore alignments and exported RDP event tables without rerunning an analysis.",className="lead text-secondary"),
    dbc.Alert(["For detection and event review, use ",html.A("nextRDP Web",href="https://murrellgroup.github.io/nextRDPweb/",target="_blank"),"."],color="info"),
    dbc.Card(dbc.CardBody([
        dcc.Upload(id="upload",children=html.Div(["Drop FASTA, CSV/TSV, or RDP5 project file here — or ",html.A("browse")]),
                   className="border rounded p-5 text-center"),
        html.Div(id="status",className="mt-3")
    ]),className="mb-3"),
    dbc.Tabs([
        dbc.Tab([dcc.Graph(id="primary-chart"),dcc.Graph(id="secondary-chart")],label="Visual overview"),
        dbc.Tab([html.Div(id="summary-cards",className="my-3"),dash_table.DataTable(id="table",page_size=15,sort_action="native",filter_action="native",
                 style_table={"overflowX":"auto"},style_cell={"fontFamily":"system-ui","fontSize":13,"padding":"7px"})],label="Data & QC"),
        dbc.Tab(dbc.Card(dbc.CardBody([
            html.H4("Supported inputs"),
            dcc.Markdown("""- **FASTA alignments:** sequence-level gap/ambiguity QC and sliding-window variability.\n- **CSV/TSV event tables:** breakpoint intervals, method support, p-values, and sortable records.\n- **RDP5 projects:** file validation and label inventory. The binary project format is not publicly documented, so event-level parsing is intentionally not claimed yet."""),
            html.H4("Privacy",className="mt-3"),html.P("Files are processed in this Dash session and are not sent to WebRDP by this app.")
        ]),className="my-3"),label="About formats")
    ]),
    html.P("Teaching and exploratory visualization only; inspect primary evidence before making biological conclusions.",className="small text-secondary text-center my-4")
],fluid="xl")


@callback(Output("parsed","data"),Output("status","children"),Input("upload","contents"),State("upload","filename"),prevent_initial_call=True)
def ingest(contents,filename):
    if not contents:return no_update,no_update
    try:
        raw=decode_upload(contents); lower=(filename or "").lower()
        if lower.endswith((".fas",".fasta",".fa",".fna")):
            names,seqs=parse_fasta(raw); qc,profile=alignment_summary(names,seqs)
            data={"kind":"alignment","filename":filename,"qc":qc.to_dict("records"),"profile":profile.to_dict("records"),"n":len(names),"length":len(seqs[0])}
            msg=f"Loaded {len(names):,} aligned sequences × {len(seqs[0]):,} sites."
        elif lower.endswith(".rdp5"):
            labels=inspect_rdp5(raw); data={"kind":"project","filename":filename,"labels":labels.to_dict("records"),"bytes":len(raw)}
            msg=f"Validated an RDP5 project ({len(raw)/1e6:.1f} MB); extracted {len(labels):,} probable labels."
        else:
            events=normalize_events(raw); data={"kind":"events","filename":filename,"events":events.to_dict("records")}
            msg=f"Loaded {len(events):,} candidate event rows."
        return data,dbc.Alert(msg,color="success",className="py-2")
    except Exception as exc:return no_update,dbc.Alert(str(exc),color="danger",className="py-2")


@callback(Output("primary-chart","figure"),Output("secondary-chart","figure"),Output("table","data"),Output("table","columns"),Output("summary-cards","children"),Input("parsed","data"))
def render(data):
    blank=px.scatter(title="Upload a lesson file or exported event table")
    if not data:return blank,blank,[],[],""
    kind=data["kind"]
    if kind=="alignment":
        qc=pd.DataFrame(data["qc"]); profile=pd.DataFrame(data["profile"])
        primary=px.line(profile,x="Midpoint",y=["Variable sites","Mean entropy"],title="Sliding-window sequence variability")
        secondary=px.bar(qc.sort_values("Gap %",ascending=False),x="Sequence",y=["Gap %","Ambiguous %"],barmode="group",title="Sequence quality overview")
        summary=dbc.Row([dbc.Col(dbc.Alert(f'{data["n"]:,} sequences',color="primary")),dbc.Col(dbc.Alert(f'{data["length"]:,} sites',color="secondary")),dbc.Col(dbc.Alert(f'{int(profile["Variable sites"].sum()):,} variable-site counts',color="info"))])
        frame=qc
    elif kind=="events":
        frame=pd.DataFrame(data["events"]); frame["Length"]=(frame["End"]-frame["Start"]).abs()+1
        primary=px.bar(frame,x="Length",y="Recombinant",base="Start",color="Method",orientation="h",hover_data=list(frame.columns),title="Candidate recombination intervals")
        counts=frame.groupby(["Recombinant","Method"]).size().reset_index(name="Events")
        secondary=px.bar(counts,x="Recombinant",y="Events",color="Method",title="Method support by recombinant")
        summary=dbc.Alert(f'{len(frame):,} events across {frame["Method"].nunique():,} methods and {frame["Recombinant"].nunique():,} recombinant labels',color="info")
    else:
        frame=pd.DataFrame(data["labels"]); primary=px.histogram(frame.assign(Name_length=frame.iloc[:,0].str.len()),x="Name_length",title="Extracted label-length distribution")
        secondary=px.scatter(title="Event plots require an exported CSV/TSV table")
        summary=dbc.Alert(f'RDP5 signature valid · {data["bytes"]/1e6:.1f} MB · {len(frame):,} probable labels',color="warning")
    for fig in (primary,secondary):fig.update_layout(template="plotly_white")
    return primary,secondary,frame.to_dict("records"),[{"name":c,"id":c} for c in frame.columns],summary


if __name__=="__main__":app.run(debug=True,host="0.0.0.0",port=8050)
