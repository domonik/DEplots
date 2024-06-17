import pandas as pd
import plotly.graph_objs as go
import plotly.express as px
from DEplots.dashboard import get_data
from plotly.subplots import make_subplots
from typing import Tuple, List
import numpy as np



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



if __name__ == '__main__':
    config_file = "/home/rabsch/PythonProjects/DEPlots/testData/config.yaml"
    rd = "/home/rabsch/PythonProjects/RlocSeq/Pipeline/RUNS/"
    _, data = get_data(config_file, rd)
    data = data[list(data.keys())[0]]
    genes = data[(data[("Additional Data", "gene_name")].str.contains("psb") == True) & (~data.index.str.contains("UTR"))].index

    runs = ["LightRunMinusPuromycin", "LightDark_LightOnly", "LightDark_DarkOnly"]
    plot_gene_among_conditions(data, genes[11:14], name_col=("Additional Data", "gene_name"), runs=None)
