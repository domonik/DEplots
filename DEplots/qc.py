from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram
import pandas as pd
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import numpy as np
import plotly.express as px


def sample_pca(df, colors = None):
    group_cols = [col for col in df.columns if col not in ["PC1", "PC2", "group", "name"]]
    if len(group_cols) == 2:
        color_column = group_cols[0]
        shape_col = group_cols[1]

    else:
        color_column = "group"
        shape_col = None
    if colors:
        unique_groups = df[color_column].unique().tolist()
        assert len(unique_groups) <= len(colors), "Not enough colors specified. there are more groups than colors"
        custom_colors = {group: colors[i] for i, group in enumerate(unique_groups)}
    else:
        custom_colors = None
    fig = px.scatter(df, y="PC2", x="PC1", color=color_column, symbol=shape_col, hover_data=["name"], color_discrete_map=custom_colors)
    return fig




def pheatmap(df, colorscale = None, tree_color: str = "black", **kwargs):
    if df.iloc[0, 0] == 1:
        matrix = 1 - df.values
    else:
        matrix = df.values
    condensed_dist = squareform(matrix, checks=False)
    Z = linkage(condensed_dist, method="complete")
    ddata = dendrogram(Z, no_plot=True)
    x_coords = (np.asarray(ddata['icoord']) / 10) - .5
    y_coords = np.asarray(ddata['dcoord'])
    default_args = {
        "row_heights": [0.2, 0.8],
        "column_widths": [0.8, 0.2],
        "shared_xaxes": "columns",
        "shared_yaxes": "rows",

    }
    kwargs = default_args | kwargs
    fig = make_subplots(rows=2, cols=2, **kwargs)
    for i in range(len(x_coords)):
        fig.add_trace(
            go.Scatter(
                mode="lines",
                x=x_coords[i],
                y=y_coords[i],
                marker=dict(color=tree_color),
                line=dict(color=tree_color),
                showlegend=False
            ),
            row=1, col=1

        )
        fig.add_trace(
            go.Scatter(
                mode="lines",
                x=y_coords[i],
                y=x_coords[i],
                marker=dict(color=tree_color),
                line=dict(color=tree_color),
                showlegend=False

            ),
            row=2, col=2

        )
    c_order = ddata["leaves"]
    heatmap_matrix = (df.iloc[c_order, :].iloc[:, c_order]).values
    ticks = df.columns[c_order]
    fig.add_trace(
        go.Heatmap(
            z=heatmap_matrix,
            x=ticks,
            y=ticks,
            colorscale=colorscale,

        ),
        row=2, col=1
    )
    fig.update_xaxes(
        range=[-.5, len(matrix) -.5],
        row=1, col=1
    )
    fig.update_yaxes(
        range=[-.5, len(matrix) - .5],
        col=2, row=2
    )
    fig.update_yaxes(
        range=[0, y_coords.max() * 1.05],
        row=1
    )
    fig.update_xaxes(
        range=[0, y_coords.max() * 1.05],
        col=2
    )
    return fig



if __name__ == '__main__':
    file = "/home/rabsch/PythonProjects/RlocSeq/Pipeline/RUNS/LightRunMinusPuromycin/PipelineData/IntermediateData/CorrData.tsv"
    file = "/home/rabsch/PythonProjects/RlocSeq/Pipeline/RUNS/LightDark_LightTotalCell/PipelineData/IntermediateData/PCAData.tsv"
    df = pd.read_csv(file, sep="\t")
    sample_pca(df, colors= ["orange", "blue"])
    pheatmap(df, colorscale= ["orange", "blue"])
