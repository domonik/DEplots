import plotly.graph_objs as go
import pandas as pd
import numpy as np
from typing import Dict, Tuple

HOVERTEMPLATE = '<i>Y</i>: %{y:.2f}' + \
                '<br><b>X</b>: %{x}<br>' + \
                '<b>%{text}</b>'

def _internal_plot_fct(df, x: str = "log2FoldChange", y: str = "-log10padj", highlight = None, name_col: str = None):
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
                x=to_highlight[x],
                y=to_highlight[y],
                mode="markers",
                marker=dict(color=color),
                hovertemplate=HOVERTEMPLATE,
                text=to_highlight[name_col] if name_col else to_highlight.index,
                name=key,
                showlegend=True
            ))
    if len(not_highlighted) >= 1:
        remaining = df[np.logical_and.reduce(not_highlighted)]
    else:
        remaining = df
    fig.add_trace(go.Scatter(
        x=remaining[x],
        y=remaining[y],
        mode="markers",
        marker=dict(color="grey"),
        hovertemplate=HOVERTEMPLATE,
        text=remaining[name_col] if name_col else remaining.index,
        name="non Highlighted Genes",
        showlegend=True
    ))
    fig.data = fig.data[::-1]
    return fig

def ma_from_deseq_result(
        deseq_result: pd.DataFrame,
        name_col: str = None,
        highlight: Dict[str, Tuple[str, str]] = None,
        lfc_cutoff: float = None,
        padj_cutoff: float = None,
        condition_name: str = None,
        base_name: str = None,
        highlight_up_color: str = "#2ca02c",
        highlight_down_color: str = "#002695",
):
    df = deseq_result[~pd.isna(deseq_result["baseMean"])]
    df["log10BaseMean"] = np.log10(df["baseMean"])
    if padj_cutoff is not None:
        lfc_cutoff = lfc_cutoff if lfc_cutoff else 0
        upreg = df[(df["padj"] <= padj_cutoff) & (df["log2FoldChange"] > lfc_cutoff)]
        downreg = df[(df["padj"] <= padj_cutoff) & (df["log2FoldChange"] < -lfc_cutoff)]
        condition_name = condition_name if condition_name else "Up"
        base_name = base_name if base_name else "Down"
        reg = {condition_name: [highlight_up_color, upreg.index.tolist()], base_name: [highlight_down_color, downreg.index.tolist()]}
        highlight = (highlight | reg) if highlight else reg
    fig = _internal_plot_fct(df, x="log10BaseMean", y="log2FoldChange", name_col=name_col, highlight=highlight)

    fig.update_layout(
        yaxis_title="Log<sub>2</sub>(FoldChange)",
        xaxis_title="Log<sub>10</sub>(mean of normalized counts)"
    )
    return fig




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

    df = deseq_result[~pd.isna(deseq_result["padj"])]
    df["-log10padj"] = -1 * np.log10(df["padj"])
    max_log10padj = np.ceil(df["-log10padj"].max())
    min_fc = np.floor(df["log2FoldChange"].min())
    max_fc = np.ceil(df["log2FoldChange"].max())

    fig = _internal_plot_fct(df, name_col=name_col, highlight=highlight)

    if lfc_cutoff is not None and padj_cutoff is not None:
        fig = add_boxes(fig, lfc_cutoff, padj_cutoff, highlight_up_color, highlight_down_color, opacity)
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
        xaxis_title="Log<sub>2</sub>(FoldChange)",
        yaxis_title="-Log<sub>10</sub>(pval)"
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