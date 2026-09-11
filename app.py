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
    # RDP5 stores the sequence names near the beginning as fixed-width,
    # NUL-padded ASCII fields.  Stop at the first long nucleotide run so the
    # embedded alignment is not mistaken for thousands of project labels.
    alignment_start=re.search(rb"[ACGTN?-]{250,}",raw)
    header=raw[:alignment_start.start()] if alignment_start else raw[:2_000_000]
    # Sequence names are stored as 100-byte, space-padded fields introduced
    # by the two-byte marker ``d\0``.  Parsing the fields directly also keeps
    # perfectly valid short lesson labels such as A, B, ... Y.
    fields=re.findall(rb"d\x00([^\x00]{1,100})",header)
    seen=[]
    for field in fields:
        value=field.decode("ascii",errors="ignore").strip()
        if value and re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_. -]{0,99}",value) and value not in seen:
            seen.append(value)
    # Fallback for project variants that do not use the fixed-width marker.
    if not seen:
        chunks=re.findall(rb"[A-Za-z0-9][A-Za-z0-9_.-]{5,95}",header)
        for chunk in chunks:
            value=chunk.decode("ascii",errors="ignore").strip("._-")
            if ("." in value or "_" in value) and value not in seen and not value.isdigit(): seen.append(value)
    rows=[]
    for label in seen[:5000]:
        parts=label.split(".")
        first=parts[0] if parts else ""
        synthetic=re.match(r"SN\d+_([A-Za-z0-9]+)$",first,re.I)
        subtype=synthetic.group(1) if synthetic else ("Unknown" if first.upper().startswith("SN") else first)
        country=parts[1] if len(parts)>2 and 2 <= len(parts[1]) <= 3 else "Unknown"
        year=parts[2] if len(parts)>3 and re.fullmatch(r"(?:\d{2}|19\d{2}|20\d{2}|x)",parts[2],re.I) else "Unknown"
        if re.fullmatch(r"\d{2}",year):
            value=int(year); year=str(1900+value if value>=70 else 2000+value)
        elif year.lower()=="x": year="Unknown"
        rows.append({"Sequence label":label,"Subtype":subtype or "Unknown","Country code":country,"Year":year})
    columns=["Sequence label","Subtype","Country code","Year"]
    return pd.DataFrame(rows,columns=columns).drop_duplicates("Sequence label")


def rdp5_alignment_overview(raw: bytes, labels: pd.DataFrame, bins: int = 120, max_tracks: int = 72) -> dict:
    """Build an RDP-like, browser-sized overview from the embedded alignment."""
    runs=[m.group().decode("ascii") for m in re.finditer(rb"[ACGTN?.-]{500,}",raw)]
    if not runs:
        return {"matrix":[],"names":[],"x":[],"diversity":[],"gaps":[],"length":0,"sequences":0}
    length=min(map(len,runs)); arr=np.array([list(s[:length]) for s in runs],dtype="U1")
    called=np.isin(arr,list("ACGT"))
    consensus=[]
    for col in arr.T:
        counts=Counter(col[np.isin(col,list("ACGT"))]); consensus.append(counts.most_common(1)[0][0] if counts else "N")
    consensus=np.asarray(consensus)
    edges=np.linspace(0,length,bins+1,dtype=int)
    matrix=np.zeros((len(runs),bins)); diversity=[]; gaps=[]
    for j,(left,right) in enumerate(zip(edges[:-1],edges[1:])):
        block=arr[:,left:right]; valid=called[:,left:right]
        mismatch=(block!=consensus[left:right]) & valid
        matrix[:,j]=100*mismatch.sum(axis=1)/np.maximum(valid.sum(axis=1),1)
        diversity.append(float(np.mean(matrix[:,j])))
        gaps.append(float(100*np.mean(~valid)))
    order=np.argsort(matrix.mean(axis=1))[::-1]
    keep=order[np.linspace(0,len(order)-1,min(max_tracks,len(order)),dtype=int)]
    label_names=labels.get("Sequence label",pd.Series(dtype=str)).tolist()
    names=[label_names[i] if i<len(label_names) else f"Sequence {i+1}" for i in keep]
    return {"matrix":np.round(matrix[keep],2).tolist(),"names":names,
            "x":[int((a+b)/2)+1 for a,b in zip(edges[:-1],edges[1:])],
            "diversity":np.round(diversity,2).tolist(),"gaps":np.round(gaps,2).tolist(),
            "length":length,"sequences":len(runs)}


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
        dbc.Tab([dcc.Graph(id="primary-chart",config={"displaylogo":False}),dcc.Graph(id="secondary-chart",config={"displaylogo":False})],label="Visual overview"),
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
            labels=inspect_rdp5(raw); overview=rdp5_alignment_overview(raw,labels)
            data={"kind":"project","filename":filename,"labels":labels.to_dict("records"),"overview":overview,"bytes":len(raw)}
            msg=(f"Loaded RDP5 project: {overview['sequences']:,} aligned sequences × "
                 f"{overview['length']:,} sites; extracted {len(labels):,} project labels.")
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
        frame=pd.DataFrame(data["labels"])
        overview=data.get("overview",{}); x=overview.get("x",[])
        primary=go.Figure(go.Heatmap(z=overview.get("matrix",[]),x=x,y=overview.get("names",[]),
            colorscale=[[0,"#f7fbff"],[.15,"#c6dbef"],[.35,"#6baed6"],[.65,"#fdae61"],[1,"#b2182b"]],
            zmin=0,zmax=35,colorbar=dict(title="Divergence<br>from consensus (%)"),hovertemplate="%{y}<br>Site %{x:,}<br>Divergence %{z:.1f}%<extra></extra>"))
        primary.update_layout(title="RDP-style alignment overview",xaxis_title="Alignment position (nt)",yaxis_title="Representative sequence tracks",height=760)
        secondary=go.Figure()
        secondary.add_trace(go.Scatter(x=x,y=overview.get("diversity",[]),name="Mean divergence",line=dict(color="#167d9a",width=3),fill="tozeroy",fillcolor="rgba(22,125,154,.12)"))
        secondary.add_trace(go.Scatter(x=x,y=overview.get("gaps",[]),name="Uncalled / gap sites",line=dict(color="#d97706",width=2)))
        secondary.update_layout(title="Genome-wide signal profile",xaxis_title="Alignment position (nt)",yaxis_title="Percent",hovermode="x unified")
        summary=dbc.Row([
            dbc.Col(dbc.Alert(f'{overview.get("sequences",0):,} aligned sequences',color="primary")),
            dbc.Col(dbc.Alert(f'{overview.get("length",0):,} nucleotide sites',color="secondary")),
            dbc.Col(dbc.Alert(f'{frame["Subtype"].replace("Unknown",np.nan).nunique():,} subtypes',color="info")),
            dbc.Col(dbc.Alert(f'{frame["Country code"].replace("Unknown",np.nan).nunique():,} country codes',color="success"))])
    for fig in (primary,secondary):
        fig.update_layout(template="plotly_white",margin=dict(l=45,r=25,t=65,b=45),hoverlabel=dict(namelength=-1))
    return primary,secondary,frame.to_dict("records"),[{"name":c,"id":c} for c in frame.columns],summary


if __name__=="__main__":app.run(debug=True,host="0.0.0.0",port=8050)
