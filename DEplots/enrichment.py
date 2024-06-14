import plotly.express as px
import plotly.graph_objs as go
import pandas as pd
from plotly.subplots import make_subplots
import numpy as np



def empty_figure(annotation: str = None):
    fig = go.Figure()
    fig.update_yaxes(showticklabels=False, showgrid=False)
    fig.update_xaxes(showgrid=False, showticklabels=False)
    fig.update_layout(
        margin={"t": 0, "b": 0, "r": 50},
        font=dict(
            size=16,
        ),
        yaxis=dict(zeroline=False),
        xaxis=dict(zeroline=False),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",

    )
    if annotation is not None:
        fig.add_annotation(
            xref="paper",
            yref="paper",
            xanchor="center",
            yanchor="middle",
            x=0.5,
            y=0.5,
            text=annotation,
            showarrow=False,
            font=(dict(size=28))
        )
    fig.layout.template = "plotly_white"
    fig.update_layout(
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),

    )
    return fig



def enrichment_plot_from_cp_table(df, mode="scatter", colorscale = None):
    def df_div(l):
        return int(l[0]) / int(l[1])
    df["Gene Ratio"] = df["GeneRatio"].str.split("/").apply(df_div)
    if len(df) == 0:
        return empty_figure("Nothing enriched")
    if mode == "scatter":
        fig = px.scatter(
            df,
            x="Gene Ratio",
            y="Description",
            symbol="ONTOLOGY" if "ONTOLOGY" in df else None,
            color="p.adjust",
            template="plotly_white",
            color_continuous_scale=colorscale

        )
        fig.update_traces(marker=dict(size=15))

    elif mode == "bar":
        fig = px.bar(
            df,
            x="Gene Ratio",
            y="Description",
            color="p.adjust",
            template="plotly_white",
        )
    else:
        raise ValueError(f"mode: {mode} is not valid")
    categories = df["Description"].unique()
    y_range = [-0.75, len(categories)]

    # Update the layout to reduce the extra space
    fig.update_layout(
        yaxis=dict(range=y_range)
    )
    fig.update_layout(
        coloraxis_colorbar=dict(
            yanchor="top",
            y=0.7,
            len=0.7,
            x=1,
            ticks="outside"
        ),
        legend=dict(x=1),
        yaxis=dict(tickmode="linear", type="category", dtick=1)
    )
    return fig


def plot_gsea(
        df: pd.DataFrame,
        info: pd.DataFrame = None,
        descs=None, colors=px.colors.DEFAULT_PLOTLY_COLORS,
        show_zero_lfc: bool = False,
        condition_name: str = None,
        base_name: str = None,
        gene_list_name: str = None

):
    if gene_list_name:
        df = df.rename({"geneList": gene_list_name}, axis=1)
    else:
        gene_list_name = "geneList"
    if info is not None:
        df = df.merge(info, left_on="Description", right_on="Description")
        hover_data = ["p.adjust", gene_list_name]
    else:
        hover_data = [gene_list_name]
    if descs is None:
        descs = df["Description"].unique()
    d = descs
    df = df[df["Description"].isin(d)]

    df['Description'] = pd.Categorical(df['Description'], categories=d, ordered=True, )
    df = df.sort_values(by=["Description", "x"])
    fig_old = px.line(df, x="x", y="runningScore", color='Description', hover_data=hover_data)
    hovertemplate = '<i>Y</i>: %{y:.2f}' + \
                    '<br><b>X</b>: %{x}<br>' + \
                    '<b>%{text}</b>'
    to_display = min(len(d), len(colors))
    fig = make_subplots(
        rows=to_display + 1, shared_xaxes=True, vertical_spacing=0,
        row_heights=[0.6] + [0.4 / to_display for _ in range(to_display)]
    )
    df[df["position"] == 0]["position"] = np.nan
    check = df[df[gene_list_name] == 0]
    if len(check) > 0:
        pos = check["x"].mean()
    else:
        pos = (df[df[gene_list_name] > 0]["x"].max() * 2 + 1) / 2

    for idx, data in enumerate(fig_old.data):
        if idx < len(colors):
            data.update(marker=dict(color=colors[idx]), line=dict(color=colors[idx]))
            fig.add_trace(data, row=1, col=1)
    for idx, desc in enumerate(d, 1):
        if idx < len(colors):
            sdf = df[df["Description"] == desc]
            sdf = sdf[sdf["position"] != 0]
            y = sdf["position"]
            fig.add_trace(
                go.Bar(
                    x=sdf["x"],
                    y=y,
                    marker_color=colors[idx-1],
                    marker_line=dict(width=1, color=colors[idx-1]),
                    name=desc,
                    hovertext=df["gene"],
                    hovertemplate=hovertemplate,
                ),
                row=idx + 1, col=1
            )
            if show_zero_lfc:
                fig.add_vline(x=pos, line_dash="dot", row=idx+1)

            fig.update_yaxes(range=[0, 1], row=idx+1, showticklabels=False)

    if show_zero_lfc:

        fig.add_vline(
            x=pos,
            annotation=dict(text="0 L2FC"), row=1,
            line_dash="dot"
        )
    if condition_name:
        fig.add_annotation(
            text=condition_name,
            showarrow=False,
            xref="x domain",
            xanchor="left",
            yanchor="top",
            yref="y domain",
            x=0,
            y=1
        )
    if base_name:
        fig.add_annotation(
            text=base_name,
            showarrow=False,
            xref="x domain",
            xanchor="right",
            yanchor="top",
            yref="y domain",
            x=1,
            y=1
        )
    fig.update_xaxes(title_text="Position in ranked dataset", row=len(d)+1)
    fig.update_yaxes(title_text="Running enrichment score", row=1)
    return fig

if __name__ == '__main__':
    file = "../../RlocSeq/Pipeline/RUNS/LightDark_DarkOnly/PipelineData/Enrichment/GSEAGO_plot_data_cM_vs_bC.tsv"
    file2 = "../../RlocSeq/Pipeline/RUNS/LightDark_DarkOnly/PipelineData/Enrichment/GSEAGO_cM_vs_bC.tsv"
    df = pd.read_csv(file, sep="\t")
    idx =pd.read_csv(file2, sep="\t")
    fig = plot_gsea(df, idx, show_zero_lfc=True, condition_name="FOOO")
    fig.show()