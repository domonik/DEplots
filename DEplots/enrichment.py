import plotly.express as px
import plotly.graph_objs as go


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