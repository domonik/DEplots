import dash
import pandas as pd
from dash import Input, Output, State, html, dcc, callback, dash_table
import dash_bootstrap_components as dbc
from dash.dash_table.Format import Format, Scheme
from pandas.api.types import is_numeric_dtype
from DEplots.runComparison import plot_gene_among_conditions, upset_plot_from_deseq
from DEplots.dashboard import DEFAULT_PLOTLY_COLORS, DEFAULT_PLOTLY_COLORS_LIST, LAYOUT, DARK_LAYOUT, UP_COLOR_LIGHT, \
    UP_COLOR_DARK, DOWN_COLOR_LIGHT, DOWN_COLOR_DARK, DASH_DATA
from DEplots.enrichment import enrichment_plot_from_cp_table, empty_figure, plot_gsea



dash.register_page(__name__, path='/comparison', name="Comparison")


def all_comparisons():
    comps = []
    for key, value in DASH_DATA[0].items():
        comp = DASH_DATA[0][key]["comparisons"]
        comps += comp
    return list(set(comps))


def get_datasets_with_comp(comp):
    ds = []
    for key, value in DASH_DATA[0].items():
        if comp in DASH_DATA[0][key]["comparisons"]:
            ds.append(key)
    return ds


def datasets_card():
    comps = list(DASH_DATA[1].keys())
    dataset_card = dbc.Col(
        dbc.Card(
            dbc.Row(
                [
                    dbc.Col(html.H5("Datasets"), width=2, className="d-flex"),
                    dbc.Col(
                        dcc.Dropdown(
                            comps, comps[0],
                            id="comparison-hash-dd",
                            clearable=False

                        ),
                        width=3

                    ),
                    dbc.Col(
                        dcc.Dropdown(
                            id="dataset-compare-dd",
                            clearable=True,
                            multi=True

                        ),
                        width=6

                    ),


                ],
                className="m-2"

            ), className="shadow"
        ), className="mt-2", width=12

    )
    return dataset_card


def get_table():
    table = dbc.Col(dbc.Card(
        [
            dbc.CardHeader(
                dbc.Row(
                    [
                        dbc.Col(html.H5("Multi DESeq table"), width=6, align="center"),
                    ],
                    justify="between"
                ),

            ),
            dbc.Row(
                [

                    dbc.Col(
                        dash_table.DataTable(
                            id='datasets-table',
                            filter_action="native",
                            filter_options={"case": "insensitive"},
                            sort_action="native",
                            sort_mode="multi",
                            column_selectable="single",
                            row_selectable="multi",
                            row_deletable=False,
                            merge_duplicate_headers=True,
                            selected_columns=[],
                            selected_rows=[],
                            page_action="native",
                            page_current=0,
                            page_size=10,
                            style_as_list_view=True,
                            fixed_columns={'headers': True, 'data': 1},
                            style_table={'minWidth': '100%'},
                            style_data_conditional=[
                                {
                                    'if': {'row_index': 'odd'},
                                    'backgroundColor': 'var(--bs-secondary-bg)',
                                },
                                {
                                    'if': {'row_index': 'even'},
                                    'backgroundColor': 'var(--bs-tertiary-bg)',
                                },

                            ],
                            style_header={
                                'backgroundColor': 'var(--bs-secondary-bg)',
                                'fontWeight': 'bold',
                                "border": "none"
                            },
                            style_filter={
                                'backgroundColor': 'var(--bs-secondary-bg)',
                                'fontWeight': 'bold',
                                "border": "none !important"

                            },
                            style_data={'border': 'none !important'}
                        ),
                        width=12, style={"overflow": "auto", 'backgroundColor': 'var(--bs-primary-bg)'},
                    )

                ],
                justify="center", className="m-2"

            )

        ],
        className="shadow"
    ), width=12)
    return table


def get_line_plot_card():
    gsea_box = dbc.Col(
        dbc.Card(
            [
                dbc.CardHeader(
                    dbc.Row(
                        [
                            dbc.Col(html.H5("Single Gene Behaviour"), width=6),
                            dbc.Col(html.Span("Name Column"), width=3,
                                    className="d-flex align-items-center justify-content-end"),
                            dbc.Col(dcc.Dropdown(
                                style={"width": "100%"},
                                id="plot-hover-name-dd",
                                clearable=False

                            ), width=3, className="d-flex align-items-center"),
                        ]

                    ),

                ),
                dbc.Col(dcc.Graph(id="gene-line-graph", ), width=12,

                        ),

            ],
            className="shadow",
        ),
        width=12
    )
    return gsea_box


def get_upset_card(updown):
    gsea_box = dbc.Col(
        dbc.Card(
            [
                dbc.CardHeader(
                    dbc.Row(
                        [
                            dbc.Col(html.H5(f"Upset {updown} regulated"), width=6),

                        ]

                    ),

                ),
                dbc.Col(dcc.Graph(id=f"upset-{updown}-graph", ), width=12,

                        ),

            ],
            className="shadow",
        ),
        width=12
    )
    return gsea_box



def get_layout(dash_data):
    lout = html.Div(
        [
            dbc.Container(
                [
                    dbc.Row(
                        datasets_card(),
                        className="py-1"

                    ),
                    dbc.Row(
                        get_line_plot_card(),
                        className="py-1"
                    ),
                    dbc.Row(
                        get_table(),
                        className="py-1"

                    ),
                    dbc.Row(
                        [
                            get_upset_card("up"),

                        ],

                        className="py-1"

                    ),
                    dbc.Row(
                        [
                            get_upset_card("down"),

                        ],

                        className="py-1"

                    ),


                ],
                fluid=True,
                className="dbc"
            )
        ]
    )
    return lout


def upset_fig(datasets, switch, comp, updown):
    dot_color = "grey" if switch else "rgb(80,80,80)"
    if datasets is None or len(datasets) == 0:
        fig = empty_figure("No dataset selected")
    else:
        idx = pd.IndexSlice
        df = DASH_DATA[1][comp]
        df = df.loc[:, idx[datasets, :]]
        if updown == "up":
            lfc_cutoff = 0.8
            barcolor = UP_COLOR_LIGHT if switch else UP_COLOR_DARK
            dot_colors = (dot_color, barcolor)

        else:
            lfc_cutoff = -0.8
            barcolor = DOWN_COLOR_LIGHT if switch else DOWN_COLOR_DARK
            dot_colors = (dot_color, barcolor)
        fig = upset_plot_from_deseq(df, padj_cutoff=0.05, lfc_cutoff=lfc_cutoff, vertical_spacing=0, bar_color=barcolor, dot_colors=dot_colors, horizontal_spacing=0, mode=updown)
    if not switch:
        fig.update_layout(DARK_LAYOUT)
        linecolor = "white"
    else:
        fig.update_layout(LAYOUT)
        linecolor = "black"
    fig.update_yaxes(showline=True, linecolor=linecolor, mirror=True)
    fig.update_xaxes(showline=True, linecolor=linecolor, mirror=True)
    fig.update_xaxes(gridcolor=dot_color, row=2, col=1)
    return fig

@callback(
    Output("upset-up-graph", "figure"),
    Input("dataset-compare-dd", "value"),
    Input("mode-switch", "value"),
    State("comparison-hash-dd", "value"),

)
def plot_upset_up(datasets, switch, comp):
    fig = upset_fig(datasets, switch, comp, "up")
    return fig

@callback(
    Output("upset-down-graph", "figure"),
    Input("dataset-compare-dd", "value"),
    Input("mode-switch", "value"),
    State("comparison-hash-dd", "value"),

)
def plot_upset_down(datasets, switch, comp):
    fig = upset_fig(datasets, switch, comp, "down")
    return fig


@callback(
    Output("plot-hover-name-dd", "options"),
    Output("plot-hover-name-dd", "value"),
    Input("comparison-hash-dd", "value"),

)
def update_name_selection(comp):
    if comp is None:
        raise dash.exceptions.PreventUpdate
    datasets = get_datasets_with_comp(comp)
    options = [f"{col[0]} - {col[1]}" for col in DASH_DATA[1][comp].columns if col[0] not in datasets]
    return options, options[0]

@callback(
    Output("datasets-table", "columns"),
    Output("datasets-table", "data"),
    Input("dataset-compare-dd", "value"),
    State("comparison-hash-dd", "value")
)
def update_datasets_table(datasets, comp):
    if datasets is None:
        return None, None
    idx = pd.IndexSlice
    df = DASH_DATA[1][comp]
    df = df.loc[:, idx[["Name", "Additional Data"] + datasets, :]]

    columns = [
        {"name": i, "id": "_".join(i), "deletable": False, "selectable": False, "format": Format(precision=4), "type": "numeric" if is_numeric_dtype(df[i]) else "text"
         } for i in df.columns
    ]
    data = [{"_".join(col): val for col, val in row.items()} for row in df.to_dict('records')]
    return columns, data


@callback(
    Output("dataset-compare-dd", "options"),
    Output("dataset-compare-dd", "value"),
    Input("comparison-hash-dd", "value"),

)
def update_selectable_datasets(comp):
    sets = get_datasets_with_comp(comp)
    return sets, sets



@callback(
    Output("gene-line-graph", "figure"),
    Input("datasets-table", "selected_rows"),
    Input("dataset-compare-dd", "value"),
    Input("mode-switch", "value"),
    Input("plot-hover-name-dd", "value"),

    State("comparison-hash-dd", "value"),

)
def update_line_plot(sel_rows, datasets, switch, legend_name, comp):
    plot = False
    if datasets is None or len(datasets) == 0:
        fig = empty_figure("No dataset selected")
    elif sel_rows is None or len(sel_rows) == 0:
        fig = empty_figure("No gene selected")
    else:
        idx = pd.IndexSlice
        df = DASH_DATA[1][comp]
        df = df.loc[:, idx[["Name", "Additional Data"] + datasets, :]]
        genes = df.iloc[sel_rows]
        legend_name = tuple(legend_name.split(" - "))
        fig = plot_gene_among_conditions(
            genes,
            genes=genes.index,
            runs=datasets,
            name_col=legend_name,
            colors=DEFAULT_PLOTLY_COLORS_LIST
        )
        fig.update_xaxes()
        plot = True
    if not switch:
        fig.update_layout(DARK_LAYOUT)
        linecolor = "white"
    else:
        fig.update_layout(LAYOUT)
        linecolor = "black"
    if plot:
        fig.update_shapes(line=dict(color=linecolor))
        fig.update_traces(
            error_y=dict(color=linecolor)
        )

        fig.update_xaxes(showgrid=False, zeroline=False, showline=True, layer="above traces", linecolor=linecolor, mirror=True)
        fig.update_yaxes(showgrid=True, zeroline=False, showline=True, mirror=True, linecolor=linecolor, nticks=3)

    return fig


layout = get_layout(DASH_DATA[0])