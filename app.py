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
from plotly.subplots import make_subplots
from dash import Dash, Input, Output, State, callback, dash_table, dcc, html, no_update


EVENT_ALIASES = {
    "method": ["method", "program"], "start": ["start", "begin", "breakpoint_start"],
    "end": ["end", "stop", "breakpoint_end"], "recombinant": ["recombinant", "recseq", "sequence"],
    "parent1": ["parent1", "major_parent", "parent_1"], "parent2": ["parent2", "minor_parent", "parent_2"],
    "pvalue": ["pvalue", "p-value", "p_value", "probability"],
}

# Okabe-Ito-inspired colors: distinct on projectors and for common forms of
# color-vision deficiency. A parent keeps the same color in every panel.
PARENT_COLORS = ["#0072B2", "#E69F00", "#009E73", "#CC79A7"]


def empty_figure(title: str) -> go.Figure:
    """Return a stable placeholder without invoking Plotly Express internals."""
    figure=go.Figure()
    figure.update_layout(
        template="plotly_white",title=title,height=420,
        xaxis=dict(visible=False),yaxis=dict(visible=False),
        annotations=[dict(text=title,x=.5,y=.5,xref="paper",yref="paper",showarrow=False,
                          font=dict(size=17,color="#637083"))],
        margin=dict(l=35,r=35,t=65,b=35))
    return figure


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
    # The first 32-bit integer after the signature is RDP5's zero-based final
    # sequence index (24 means 25 sequences).  Later header sections may repeat
    # or add names for saved analysis objects, so keep only alignment labels.
    stored_count=int.from_bytes(raw[17:21],"little",signed=False)+1 if len(raw)>=21 else 0
    if 0 < stored_count <= len(fields):
        fields=fields[:stored_count]
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
    return {"matrix":np.round(matrix[keep],2).tolist(),"names":names,"sequences_raw":runs,
            "x":[int((a+b)/2)+1 for a,b in zip(edges[:-1],edges[1:])],
            "diversity":np.round(diversity,2).tolist(),"gaps":np.round(gaps,2).tolist(),
            "length":length,"sequences":len(runs)}


def sliding_pairwise_identity(query: str, parent: str, window: int = 300, step: int = 50) -> tuple[list[int],list[float]]:
    """Pairwise identity using only sites called A/C/G/T in both sequences."""
    length=min(len(query),len(parent)); window=max(20,min(window,length)); step=max(1,step)
    centers=[]; identity=[]
    for left in range(0,max(1,length-window+1),step):
        right=min(left+window,length); q=np.asarray(list(query[left:right])); p=np.asarray(list(parent[left:right]))
        valid=np.isin(q,list("ACGT")) & np.isin(p,list("ACGT"))
        centers.append(left+(right-left)//2+1)
        identity.append(round(100*float(np.mean(q[valid]==p[valid])),2) if valid.any() else np.nan)
    return centers,identity


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
        dbc.Tab([
            dbc.Card(dbc.CardBody([
                dbc.Row([
                    dbc.Col([dbc.Label("Query / suspected recombinant"),dcc.Dropdown(id="query-sequence",clearable=False)],md=5),
                    dbc.Col([dbc.Label("Candidate parents (up to 4)"),dcc.Dropdown(id="parent-sequences",multi=True)],md=7)
                ],className="g-3"),
                dbc.Row([
                    dbc.Col([dbc.Label("Window (nt)"),dcc.Slider(id="similarity-window",min=100,max=800,step=50,value=300,marks={100:"100",300:"300",500:"500",800:"800"})],md=8),
                    dbc.Col([dbc.Label("Step (nt)"),dcc.Dropdown(id="similarity-step",options=[25,50,100,200],value=50,clearable=False)],md=4)
                ],className="g-3 mt-1")
            ]),className="my-3"),
            dcc.Graph(id="similarity-chart",config={"displaylogo":False}),
            dcc.Graph(id="mosaic-chart",config={"displaylogo":False}),
            dbc.Alert("Similarity switches can identify regions worth investigating, but this plot is not itself a recombination test.",color="warning")
        ],label="Similarity scan"),
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
    blank=empty_figure("Upload a lesson file or exported event table")
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


@callback(Output("query-sequence","options"),Output("query-sequence","value"),
          Output("parent-sequences","options"),Output("parent-sequences","value"),Input("parsed","data"))
def similarity_controls(data):
    if not data or data.get("kind")!="project": return [],None,[],[]
    overview=data.get("overview",{}); seqs=overview.get("sequences_raw",[])
    labels=[row.get("Sequence label",f"Sequence {i+1}") for i,row in enumerate(data.get("labels",[]))]
    labels=(labels+[f"Sequence {i+1}" for i in range(len(labels),len(seqs))])[:len(seqs)]
    options=[{"label":name,"value":i} for i,name in enumerate(labels)]
    return options,(0 if options else None),options,list(range(1,min(4,len(options))))


@callback(Output("similarity-chart","figure"),Output("mosaic-chart","figure"),Input("parsed","data"),Input("query-sequence","value"),
          Input("parent-sequences","value"),Input("similarity-window","value"),Input("similarity-step","value"))
def render_similarity(data,query_index,parent_indices,window,step):
    blank=empty_figure("Upload an RDP5 project to compare sequence similarity")
    blank_mosaic=empty_figure("Select at least two candidate parents to build a mosaic map")
    if not data or data.get("kind")!="project" or query_index is None:return blank,blank_mosaic
    overview=data.get("overview",{}); seqs=overview.get("sequences_raw",[])
    labels=[row.get("Sequence label",f"Sequence {i+1}") for i,row in enumerate(data.get("labels",[]))]
    labels=(labels+[f"Sequence {i+1}" for i in range(len(labels),len(seqs))])[:len(seqs)]
    if not 0<=int(query_index)<len(seqs):return blank,blank_mosaic
    parents=[int(i) for i in (parent_indices or []) if int(i)!=int(query_index) and 0<=int(i)<len(seqs)][:4]
    fig=go.Figure(); curves=[]; positions=[]
    for color_index,i in enumerate(parents):
        x,y=sliding_pairwise_identity(seqs[int(query_index)],seqs[i],int(window or 300),int(step or 50))
        positions=x; curves.append(y)
        hover="Site %{x:,}<br>Identity %{y:.2f}%<extra>"+labels[i]+"</extra>"
        fig.add_trace(go.Scatter(x=x,y=y,mode="lines",name=labels[i],
            line=dict(width=3,color=PARENT_COLORS[color_index]),hovertemplate=hover))
    fig.update_layout(template="plotly_white",title=f"Sliding-window similarity to {labels[int(query_index)]}",
        xaxis_title="Alignment position (nt)",yaxis_title="Pairwise identity (%)",hovermode="x unified",
        legend_title="Candidate parent",height=620,margin=dict(l=55,r=25,t=70,b=50))
    scores=np.asarray(curves,dtype=float)
    finite=scores[np.isfinite(scores)] if curves else np.array([])
    lower=max(0,float(np.floor(np.nanpercentile(finite,2)-2))) if finite.size else 0
    fig.update_yaxes(range=[lower,100])
    if len(curves)<2:return fig,blank_mosaic
    safe=np.where(np.isnan(scores),-np.inf,scores)
    raw_winners=np.argmax(safe,axis=0)
    # Suppress one-window flicker caused by nearly tied, overlapping windows.
    winners=raw_winners.copy()
    for j in range(len(winners)):
        local=raw_winners[max(0,j-2):min(len(winners),j+3)]
        winners[j]=Counter(local).most_common(1)[0][0]
    ranked=np.sort(safe,axis=0)
    confidence=ranked[-1]-ranked[-2]
    changes=np.where(winners[1:]!=winners[:-1])[0]+1
    palette=PARENT_COLORS
    for change in changes:
        fig.add_vline(x=positions[change],line_width=1,line_dash="dot",line_color="#475569",opacity=.65)
    half_step=max(1,int(step or 50)//2); segments=[]; start=0
    for end in list(changes)+[len(winners)]:
        idx=int(winners[start]); left=max(1,positions[start]-half_step); right=min(overview.get("length",positions[end-1]+half_step),positions[end-1]+half_step)
        segments.append({"parent_index":idx,"left":left,"right":right,
                         "mean_identity":float(np.nanmean(scores[idx,start:end])),
                         "mean_confidence":float(np.nanmean(confidence[start:end]))})
        start=end
    mosaic=make_subplots(rows=2,cols=1,shared_xaxes=True,row_heights=[.58,.42],vertical_spacing=.13,
                         subplot_titles=("Closest-parent segments","Confidence in the closest parent"))
    legend_seen=set()
    for segment in segments:
        idx=segment["parent_index"]; parent_name=labels[parents[idx]]; show=idx not in legend_seen; legend_seen.add(idx)
        custom=[[segment["left"],segment["right"],segment["mean_identity"],segment["mean_confidence"]]]
        mosaic.add_trace(go.Bar(x=[segment["right"]-segment["left"]+1],y=["Query mosaic"],base=[segment["left"]],
            orientation="h",width=.55,name=parent_name,legendgroup=str(idx),showlegend=show,marker_color=palette[idx],customdata=custom,
            hovertemplate="Sites %{customdata[0]:,}–%{customdata[1]:,}<br>Closest parent: "+parent_name+
                          "<br>Mean identity: %{customdata[2]:.2f}%<br>Mean lead: %{customdata[3]:.2f} pp<extra></extra>"),row=1,col=1)
    mosaic.add_trace(go.Scatter(x=positions,y=confidence,mode="lines",name="Identity lead",showlegend=False,
        line=dict(color="#334155",width=2),fill="tozeroy",fillcolor="rgba(51,65,85,.14)",
        hovertemplate="Site %{x:,}<br>Best-parent lead: %{y:.2f} percentage points<extra></extra>"),row=2,col=1)
    for change in changes:
        mosaic.add_vline(x=positions[change],line_width=2,line_color="#172033")
    mosaic.add_hline(y=1,row=2,col=1,line_dash="dot",line_color="#d97706",annotation_text="1 pp lead")
    mosaic.update_layout(template="plotly_white",title="Candidate recombination mosaic · parent switches are exploratory",
        height=470,margin=dict(l=80,r=25,t=90,b=55),barmode="overlay",
        legend=dict(orientation="h",y=1.16,x=1,xanchor="right"),hovermode="x unified")
    mosaic.update_xaxes(title_text="Alignment position (nt)",row=2,col=1)
    mosaic.update_yaxes(fixedrange=True,row=1,col=1)
    mosaic.update_yaxes(title_text="Lead (pp)",rangemode="tozero",row=2,col=1)
    return fig,mosaic


if __name__=="__main__":app.run(debug=True,host="0.0.0.0",port=8050)
