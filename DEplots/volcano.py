import plotly.graph_objs as go
import pandas as pd
import numpy as np
from typing import Dict, Tuple


def volcano_from_deseq_result(
        deseq_result: pd.DataFrame,
        name_col: str = None,
        highlight: Dict[str, Tuple[str, str]] = None,
        lfc_cutoff: float = None,
        padj_cutoff: float = None,
        condition_name: str = None,
        base_name: str = None,
        highlight_up_color: str = "#2ca02c",
        highlight_down_color: str = "#002695",
        opacity: float = 0.05
):
    hovertemplate = '<i>Y</i>: %{y:.2f}' + \
                    '<br><b>X</b>: %{x}<br>' + \
                    '<b>%{text}</b>'
    df = deseq_result[~pd.isna(deseq_result["padj"])]
    df["-log10padj"] = -1 * np.log10(df["padj"])
    max_log10padj = np.ceil(df["-log10padj"].max())
    min_fc = np.floor(df["log2FoldChange"].min())
    max_fc = np.ceil(df["log2FoldChange"].max())

    fig = go.Figure()
    not_highlighted = []
    if highlight is not None:
        for key, value in highlight.items():
            color, names = value
            mask = df.index.isin(names)
            to_highlight = df[mask]
            n = ~mask
            not_highlighted.append(n)
            fig.add_trace(go.Scatter(
                x=to_highlight["log2FoldChange"],
                y=to_highlight["-log10padj"],
                mode="markers",
                marker=dict(color=color),
                hovertemplate=hovertemplate,
                text=to_highlight[name_col] if name_col else to_highlight.index,
                name=key,
                showlegend=True
            ))
    if len(not_highlighted) >= 1:
        remaining = df[np.logical_and.reduce(not_highlighted)]
    else:
        remaining = df
    fig.add_trace(go.Scatter(
        x=remaining["log2FoldChange"],
        y=remaining["-log10padj"],
        mode="markers",
        marker=dict(color="grey"),
        hovertemplate=hovertemplate,
        text=remaining[name_col] if name_col else remaining.index,
        name="non Highlighted Genes",
        showlegend=True
    ))
    fig.data = fig.data[::-1]
    if lfc_cut_off is not None and padj_cutoff is not None:
        fig = add_boxes(fig, lfc_cut_off, padj_cutoff, highlight_up_color, highlight_down_color, opacity)
    fig.update_xaxes(
        range=[min_fc, max_fc],
        dtick=1,
        autorangeoptions=dict(
            clipmax=max_fc,
            clipmin=min_fc
        )
    )
    fig.update_yaxes(
        range=[-0.5, max_log10padj],
        autorangeoptions=dict(
            clipmax=max_log10padj,
            clipmin=-0.5
        )
    )
    fig.update_layout(
        xaxis_title="Log2FoldChange",
        yaxis_title="-Log10(pval)"
    )

    if condition_name:
        fig.add_annotation(
            text=condition_name,
            font=dict(color=highlight_up_color),
            showarrow=False,
            x=1,
            y=1,
            xref="x domain",
            yref="y domain",
            xanchor="right",
            yanchor="top"
        )
    if base_name:
        fig.add_annotation(
            text=base_name,
            font=dict(color=highlight_down_color),
            showarrow=False,
            x=0,
            y=1,
            xref="x domain",
            yref="y domain",
            xanchor="left",
            yanchor="top"
        )
    return fig


def add_boxes(fig: go.Figure, lfc_cut_off: float, padj_cutoff: float, up_color: str, down_color: str, opacity):
    fig.add_shape(
        type="rect",
        x0=lfc_cut_off, y0=-1 * np.log10(padj_cutoff), x1=100, y1=200,
        line=dict(
            color="grey",
            width=2,
        ),
        layer="below",
        fillcolor=up_color,
        opacity=opacity
    )
    fig.add_shape(
        type="rect",
        x0=-lfc_cut_off, y0=-1 * np.log10(padj_cutoff), x1=-100, y1=200,
        line=dict(
            color="grey",
            width=2,
        ),
        layer="below",
        fillcolor=down_color,
        opacity=opacity
    )
    fig.add_hline(y=-1 * np.log10(padj_cutoff), line_color="grey", layer="below", line_dash="dash")
    fig.add_vline(x=-lfc_cut_off, line_color="grey", layer="below", line_dash="dash")
    fig.add_vline(x=lfc_cut_off, line_color="grey", layer="below", line_dash="dash")
    return fig