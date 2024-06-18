import pandas as pd
import plotly.graph_objs as go
import plotly.express as px
from DEplots.dashboard import get_data
from plotly.subplots import make_subplots
from typing import Tuple, List
import numpy as np
import itertools



def plot_gene_among_conditions(df, genes, name_col: Tuple = None, runs: List = None, colors = px.colors.DEFAULT_PLOTLY_COLORS, **kwargs):

    fig = make_subplots(
        rows=len(genes),
        shared_xaxes=True,
        y_title="log2FoldChange",
        **kwargs
    )
    idx = pd.IndexSlice
    names = df.index if not name_col else df.loc[:, name_col]
    iidx = idx[:, ["log2FoldChange", "lfcSE"]] if runs is None else idx[runs, ["log2FoldChange", "lfcSE"]]

    sdf = df.loc[genes, iidx]
    columns = sdf.columns.get_level_values(0)[sdf.columns.get_level_values(1) == 'log2FoldChange'].unique()
    for condition in columns:
        sdf[(condition, 'log2FoldChange+SE')] = sdf[(condition, 'log2FoldChange')] + sdf[(condition, 'lfcSE')]
        sdf[(condition, 'log2FoldChange-SE')] = sdf[(condition, 'log2FoldChange')] - sdf[(condition, 'lfcSE')]
    max_val = np.ceil(sdf.max().max())
    min_val = np.floor(sdf.min().min())
    for i, gene in enumerate(genes, 1):
        slice = sdf.loc[gene, pd.IndexSlice[:, "log2FoldChange"]]
        errors = sdf.loc[gene, pd.IndexSlice[:, "lfcSE"]]
        name = names.loc[gene] if name_col else gene
        fig.add_trace(
            go.Bar(
                x=columns,
                y=slice,
                name=name,
                error_y=dict(
                    type="data",
                    array=errors,
                    visible=True,
                ),
                showlegend=True,
                #line=dict(color=colors[i-1]),
                marker=dict(color=colors[i-1])
            ),
            row=i, col=1

        )
        fig.update_yaxes(range=[min_val, max_val], row=i)
    fig.add_hline(y=0, line=dict(dash="dot"))
    return fig


def plotly_upset_plot(df, sorted = False, bar_color="blue", dot_colors=("grey", "black"), trim_zeros: bool = True, show_sum: bool = True,  **kwargs):
    # an array of dimensions d x d*2^d possible subsets where d is the number of columns
    subsets = []
    hovertemplate = "%{hovertext}<extra></extra>"

    # the sizes of each subset (2^d array)
    subset_sizes = []
    d = len(df.columns)
    for i in range(1, d + 1):
        subsets = subsets + [list(x) for x in list(itertools.combinations(df.columns, i))]

    for s in subsets:
        curr_bool = [1] * len(df)
        for col in df.columns:
            if col in s:
                curr_bool = [x and y for x, y in zip(curr_bool, list(df.loc[:, col].copy()))]
            else:
                curr_bool = [x and not y for x, y in zip(curr_bool, list(df.loc[:, col].copy()))]
        subset_sizes.append(sum(curr_bool))
    plot_df = pd.DataFrame({'Intersection': subsets, 'Size': subset_sizes})
    plot_df["text"] = plot_df["Intersection"].apply(lambda  x: "<br>".join(x))

    if trim_zeros:
        plot_df = plot_df[plot_df["Size"] != 0]
    if sorted:
        plot_df = plot_df.sort_values(by='Size', ascending=False)
    if "column_widths" not in kwargs:
        kwargs["column_widths"] = [0.8, 0.2]
    cols = 2 if show_sum else 1
    fig = make_subplots(rows=2, cols=cols, shared_xaxes=True, shared_yaxes=True, **kwargs)
    fig.add_trace(
        go.Bar(
            x=plot_df["text"],
            y=plot_df["Size"],
            marker=dict(color=bar_color),
            name="Intersection size",
            hovertemplate="<b>%{y}<br></b>%{x}",

        )
    )

    x = np.tile(plot_df["text"], len(df.columns))
    y = np.repeat(df.columns, len(plot_df["text"]))

    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            #line=None,
            mode="markers",
            marker=dict(color=dot_colors[0], size=20),
            showlegend=False,
            hoverinfo="skip"
        ),
        row=2,
        col=1
    )
    for idx, row in plot_df.iterrows():
        y = row["Intersection"]

        fig.add_trace(
            go.Scatter(
                x=np.repeat(row["text"], len(y)),
                y=y,
                mode="lines+markers",
                marker=dict(color=dot_colors[1], size=20, line=dict(color=dot_colors[1], width=2)),
                line=dict(color=dot_colors[1], width=5),
                showlegend=False,
                hovertext=np.repeat(row["text"], len(y)),
                hovertemplate=hovertemplate,
            ),
            row=2,
            col=1
        )
    y_range = [-0.5, len(df.columns) - 0.5]
    x_range = [-1, len(plot_df["text"])]
    if show_sum:
        df_sum = df.sum(axis=0)
        fig.add_trace(
            go.Bar(
                x=df_sum,
                y=df_sum.keys(),
                marker=dict(color=bar_color),
                name="Set size",
                orientation="h",
                hovertemplate="<b>%{x}<br></b>%{y}"

            ),
            row=2,
            col=2
        )
    # Update the layout to reduce the extra space
    fig.update_yaxes(range=y_range, row=2, showgrid=False)
    fig.update_xaxes(range=x_range, col=1)
    fig.update_traces(showlegend=False)
    fig.update_xaxes(showticklabels=True, col=2)
    fig.update_xaxes(showticklabels=False, col=1)
    fig.update_layout(hoverdistance=100)
    fig.add_annotation(
        text="Set size",
        xref="x4 domain",
        yref="y4 domain",
        xanchor="center",
        yanchor="bottom",
        showarrow=False,
        x=0.5,
        y=1
    )
    fig.add_annotation(
        text="Intersection size",
        xref="x domain",
        yref="y domain",
        xanchor="left",
        yanchor="middle",
        textangle=90,
        showarrow=False,
        x=1,
        y=0.5
    )

    #fig.update_layout(hovermode="x")
    return fig



def upset_plot_from_deseq(df, padj_cutoff, lfc_cutoff, mode: str = "up", **kwargs):
    columns = df.columns.get_level_values(0)[df.columns.get_level_values(1) == 'log2FoldChange'].unique()
    data = {}
    if mode == "up":
        for condition in columns:
            data[condition] = (df[(condition, 'log2FoldChange')] >= lfc_cutoff) & (df[(condition, 'padj')] <= padj_cutoff)
    elif mode == "down":
        for condition in columns:
            data[condition] = (df[(condition, 'log2FoldChange')] <= lfc_cutoff) & (df[(condition, 'padj')] <= padj_cutoff)

    data = pd.DataFrame(data)
    fig = plotly_upset_plot(data, **kwargs)
    return fig




if __name__ == '__main__':
    config_file = "/home/rabsch/PythonProjects/DEPlots/testData/config.yaml"
    rd = "/home/rabsch/PythonProjects/RlocSeq/Pipeline/RUNS/"
    _, data = get_data(config_file, rd)
    data = data[list(data.keys())[0]]
    genes = data[(data[("Additional Data", "gene_name")].str.contains("psb") == True) & (~data.index.str.contains("UTR"))].index

    runs = ["LightRunMinusPuromycin", "LightDark_LightOnly", "LightDark_DarkOnly"]
    fig = upset_plot_from_deseq(data, padj_cutoff=0.05, lfc_cutoff=0.8, mode="down")

    fig.show()