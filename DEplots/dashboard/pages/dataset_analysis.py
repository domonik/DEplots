from DEplots.diffexp import volcano_from_deseq_result, ma_from_deseq_result
from DEplots.enrichment import enrichment_plot_from_cp_table, empty_figure, plot_gsea
from DEplots.qc import pheatmap, sample_pca
import dash
from dash import callback, html, clientside_callback, Input, Output, dcc, dash_table, State, Patch
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
from dash.dash_table.Format import Format, Scheme
from pandas.api.types import is_numeric_dtype
from DEplots.dashboard import DEFAULT_PLOTLY_COLORS, DEFAULT_PLOTLY_COLORS_LIST, LAYOUT, DARK_LAYOUT, UP_COLOR_LIGHT, \
    UP_COLOR_DARK, DOWN_COLOR_LIGHT, DOWN_COLOR_DARK, DASH_DATA
import pandas as pd

dash.register_page(__name__, path='/analysis', name="Visualization")


def get_deseq_result(dataset_key, comp):
    df = DASH_DATA[1][comp]
    idx = pd.IndexSlice
    df1_columns = df.loc[:, idx[["Name", "Additional Data", dataset_key], :]]
    df1_columns.columns = df1_columns.columns.droplevel(0)
    return df1_columns


def get_enrich_result(dataset_key, comp, enrich: str = "GO", updown: str = "up-regulated"):
    df = DASH_DATA[0][dataset_key]["comparisons"][comp]["enrich"][enrich][updown]
    return df

def get_gsea_result(dataset_key, comp):
    df = DASH_DATA[0][dataset_key]["comparisons"][comp]["gsea"]["df"]
    plot_data = DASH_DATA[0][dataset_key]["comparisons"][comp]["gsea"]["plot_data"]
    return df, plot_data




def get_table(dash_data):
    d = list(dash_data.keys())[0]
    data = dash_data[d]
    comp = list(data["comparisons"].keys())[0]
    df = get_deseq_result(d, comp)

    table = dbc.Card(
        [
            dbc.CardHeader(
                dbc.Row(
                    [
                        dbc.Col(html.H5("DESeq table"), width=6, align="center"),
                        dbc.Col([
                            dbc.Button("Select Filtered", id="select-all", className="m-1"),
                            dbc.Button("Reset Selection", id="deselect-all", className="m-1"),
                        ], width=3, align="center", className="d-flex justify-content-end"),
                    ],
                    justify="between"
                ),

            ),
            dbc.Row(
                [

                    dbc.Col(
                        dash_table.DataTable(
                            id='deseq-table',
                            columns=None,
                            data=None,
                            editable=True,
                            filter_action="native",
                            filter_options={"case": "insensitive"},
                            sort_action="native",
                            sort_mode="multi",
                            column_selectable="single",
                            merge_duplicate_headers=True,
                            row_selectable="multi",
                            row_deletable=False,
                            selected_columns=[],
                            selected_rows=[],
                            page_action="native",
                            page_current=0,
                            page_size=10,
                            style_as_list_view=True,
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



    )
    return table


def get_gsea_table(dash_data):
    d = list(dash_data.keys())[0]
    data = dash_data[d]
    comp = list(data["comparisons"].keys())[0]
    df, _ = get_gsea_result(d, comp)
    table = dbc.Card(
        [
            dbc.CardHeader(
                dbc.Row(
                    [
                        dbc.Col(html.H5("GSEA table"), width=6, align="center"),
                    ],
                    justify="between"
                ),

            ),
            dbc.Row(
                [

                    dbc.Col(
                        dash_table.DataTable(
                            id='gsea-table',
                            columns=[
                                {"name": i, "id": i, "deletable": False, "selectable": False} for i in df.columns
                            ],
                            data=df.to_dict('records'),
                            editable=True,
                            filter_action="native",
                            filter_options={"case": "insensitive"},
                            sort_action="native",
                            sort_mode="multi",
                            column_selectable="single",
                            row_selectable="multi",
                            row_deletable=False,
                            selected_columns=[],
                            selected_rows=[],
                            page_action="native",
                            page_current=0,
                            page_size=10,
                            style_as_list_view=True,
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



    )
    return table




def get_dataset_card(dash_data):
    dataset_card = dbc.Col(
        dbc.Card(
            dbc.Row(
                [
                    dbc.Col(html.H5("Dataset"), width=2, className="d-flex"),
                    dbc.Col(
                        dcc.Dropdown(
                            [
                                name for name in dash_data
                            ], list(dash_data.keys())[0],
                            id="dataset-dd",
                            clearable=False,
                            className="navbar-input"

                        ),
                        width=6

                    ),
                    dbc.Col(
                        dcc.Dropdown(
                            [
                                name for name in dash_data
                            ], list(dash_data.keys())[0],
                            id="comparison-dd",
                            clearable=False

                        ),
                        width=3

                    )

                ],
                className="m-2"

            ), className="shadow"
        ), className="mt-2", width=12

    )
    return dataset_card



def get_gsea_box():
    gsea_box = dbc.Col(
        dbc.Card(
            [
                dbc.CardHeader(
                    dbc.Row(
                        [
                            dbc.Col(html.H5("GSEA"), width=6),
                        ]

                    ),

                ),
                dbc.Col(dcc.Graph(id="gsea-graph", ), width=12,

                        ),

            ],
            className="shadow",
        ),
        width=12
    )
    return gsea_box


def _get_qc_card():
    div = [
        dbc.Col(
            dbc.Card(
                [
                    dbc.CardHeader(
                        dbc.Row(
                            [
                                dbc.Col(html.H5("Sample Correlation"), width=6),
                            ]

                        ),

                    ),
                    dbc.Col(
                        dcc.Graph(id="correlation-graph", ), width=12,

                    ),

                ],
                className="shadow",
            ),
            width=12, md=6
        ),
        dbc.Col(
            [
                dbc.Card(
                    [
                        dbc.CardHeader(
                            dbc.Row(
                                [
                                    dbc.Col(html.H5("Sample PCA"), width=6),
                                ]

                            ),

                        ),
                        dbc.Col(
                            dcc.Graph(id="pca-graph", ), width=12,
                                ),

                    ],
                    className="shadow",
                ),

            ],
            width=12, md=6
        )
    ]
    return div


def get_layout(dash_data):
    layout = html.Div([
        dcc.Store(data={}, id="volcano-highlight-ids"),
        dcc.Store(id="enrich-term"),
        dbc.Container(
            [
                dbc.Row(
                    get_dataset_card(dash_data),
                    className="py-1"
                ),
                dbc.Row(
                    [
                        dbc.Col(
                            dbc.Card(
                                [
                                    dbc.CardHeader(
                                        dbc.Row(
                                            [
                                                dbc.Col(html.H5("DEseq"), width=2),
                                                dbc.Col(dcc.Dropdown(
                                                    value="Volcano",
                                                    options=["Volcano", "MA"],
                                                    style={"width": "100%", "font-size": "1.25rem"},
                                                    id="volcano-type-dd",
                                                    clearable=False,

                                                ), width=4, className="d-flex align-items-left ",
                                                ),
                                                dbc.Col(html.Span("Name Column"), width=3,
                                                        className="d-flex align-items-center justify-content-end"),
                                                dbc.Col(dcc.Dropdown(
                                                    value="Name",
                                                    options=["Name"],
                                                    style={"width": "100%"},
                                                    id="add-name-dd",
                                                    clearable=False

                                                ), width=3, className="d-flex align-items-center"),

                                            ]

                                        ),

                                    ),
                                    dbc.Col(dcc.Graph(id="volcano-graph", ), width=12,

                                            ),

                                ],
                                className="shadow",
                            ),
                            width=12, lg=6, className="py-1"
                        ),
                        dbc.Col(
                            [
                                dbc.Card(
                                    [
                                        dbc.CardHeader(
                                            dbc.Row(
                                                [
                                                    dbc.Col(dcc.Dropdown(
                                                        value="GO",
                                                        options=["GO", "KEGG"],
                                                        style={"width": "100%", "font-size": "1.25rem"},
                                                        id="enrich-type-dd",
                                                        clearable=False,


                                                    ), width=2, className="d-flex align-items-center ",
                                                    ),
                                                    dbc.Col(html.H5("Enrichment"), width=3,
                                                            className="d-flex align-items-center justify-content-center"),
                                                    dbc.Col(dcc.Dropdown(
                                                        value="up-regulated",
                                                        options=["up-regulated", "down-regulated"],
                                                        style={"width": "100%"},
                                                        id="enrich-updown-dd",
                                                        clearable=False

                                                    ), width=3, className="d-flex align-items-center",
                                                        style={"font-size": "1.25rem"}),
                                                    dbc.Col(
                                                        html.Span(
                                                            [
                                                                html.I(className="fas fa-xl fa-question-circle fa px-2",
                                                                       id="filter-tip"),
                                                                dbc.Tooltip(
                                                                    "Clicking on the marker of a term will highlight the"
                                                                    " genes belonging to this term in the Volcano or MA plot",
                                                                    target="filter-tip"
                                                                ),
                                                            ],

                                                        ),
                                                        width=1,
                                                        className="d-flex align-items-center justify-content-end"

                                                    ),
                                                    dbc.Col(
                                                        dbc.Button(
                                                            "Reset Selection",
                                                            disabled=True,
                                                            id="reset-enrich-selection"

                                                        ),
                                                        width=3,
                                                        className="d-flex align-items-center justify-content-end"
                                                    ),

                                                ]

                                            ),

                                        ),

                                        dbc.Col(
                                            dcc.Graph(id="enrichment-graph", ), width=12,

                                        ),

                                    ],
                                    className="shadow",

                                ),

                            ],
                            width=12, lg=6, className="py-1"
                        )
                    ],
                    className="py-1"
                ),
                dbc.Row(
                    [
                        dbc.Col(
                            get_table(dash_data),
                            width=12,
                        ),

                    ],

                    className="py-1",
                ),
                dbc.Row(
                    [
                            get_gsea_box(),

                    ],

                    className="py-1",
                ),
                dbc.Row(
                    [
                        dbc.Col(
                            get_gsea_table(dash_data),
                            width=12,
                        ),

                    ],

                    className="py-1",
                ),
                dbc.Row(
                    _get_qc_card(),
                    className="py-1",
                ),

            ],
            fluid=True,
            className="dbc",
        ),

    ])
    return layout


layout = get_layout(DASH_DATA[0])



@callback(
    [Output('deseq-table', 'selected_rows')],
    [
        Input('select-all', 'n_clicks'),
        Input('deselect-all', 'n_clicks')
    ],
    [
        State('deseq-table', 'data'),
        State('deseq-table', 'derived_virtual_data'),
        State('deseq-table', 'derived_virtual_selected_rows')
    ]
)
def select_all(select_n_clicks, deselect_n_clicks, original_rows, filtered_rows, selected_rows):
    ctx = dash.callback_context.triggered[0]
    ctx_caller = ctx['prop_id']
    if filtered_rows is not None:
        if ctx_caller == 'select-all.n_clicks':
            selected_ids = [row for row in filtered_rows]
            return [[i for i, row in enumerate(original_rows) if row in selected_ids]]
        if ctx_caller == 'deselect-all.n_clicks':
            return [[]]
        raise PreventUpdate
    else:
        raise PreventUpdate

@callback(
    Output('gsea-graph', 'figure'),
    Input('gsea-table', 'selected_rows'),
    Input('mode-switch', 'value'),

    State('dataset-dd', 'value'),
    State('comparison-dd', 'value'),
)
def get_gsea_plot(selected_rows, switch, dataset_key, comp):
    df, plot_data = get_gsea_result(dataset_key, comp)
    descs = df["Description"].iloc[selected_rows]
    cond = DASH_DATA[0][dataset_key]["comparisons"][comp]["condition"]
    baseline = DASH_DATA[0][dataset_key]["comparisons"][comp]["baseline"]
    if df is not None and len(df) > 0:
        fig = plot_gsea(plot_data, descs=descs, colors=DEFAULT_PLOTLY_COLORS_LIST, show_zero_lfc=True, condition_name=cond, base_name=baseline, gene_list_name="log2FC", vertical_spacing=0)
    elif df is None:
        fig = empty_figure("No GSEA file found for dataset")
    else:
        fig = empty_figure("Nothing enriched in this dataset")

    if not switch:
        fig.update_layout(DARK_LAYOUT)
        linecolor = "white"
    else:
        fig.update_layout(LAYOUT)
        linecolor = "black"
    fig.update_shapes(line=dict(color=linecolor))



    fig.update_xaxes(showgrid=False, zeroline=False, showline=True, layer="above traces", linecolor=linecolor)
    fig.update_yaxes(showgrid=False, zeroline=False, showline=True, mirror=True, linecolor=linecolor)

    return fig




@callback(
    Output("deseq-table", "data"),
    Output("deseq-table", "columns"),
    Output('volcano-highlight-ids', 'data'),
    Input('dataset-dd', 'value'),
    Input('comparison-dd', 'value'),
)
def update_table_from_dataset(dataset_key, comp):
    df = get_deseq_result(dataset_key, comp)
    columns = []
    for i in df.columns:
        numeric = is_numeric_dtype(df[i])
        if numeric:
            dtype = "numeric"
            if i in ("pvalue", "padj"):
                format = Format(precision=2, scheme=Scheme.exponent)
            else:
                format = Format(precision=4)
        else:
            dtype = "text"
            format = None
        d = dict(
            name=i,
            id=i,
            deletable=False,
            selectable=False,
            type=dtype,
            format=format
        )
        columns.append(d)


    data = df.to_dict('records')
    return data, columns, {}


@callback(
    Output("gsea-table", "data"),
    Output("gsea-table", "columns"),
    Output('gsea-table', 'selected_rows'),
    Input('dataset-dd', 'value'),
    Input('comparison-dd', 'value'),
)
def update_table_from_dataset(dataset_key, comp):
    df, _ = get_gsea_result(dataset_key, comp)
    columns = []
    for i in df.columns:
        numeric = is_numeric_dtype(df[i])
        if numeric:
            dtype = "numeric"
            if i in ("pvalue", "padj"):
                format = Format(precision=2, scheme=Scheme.exponent)
            else:
                format = Format(precision=4)
        else:
            dtype = "text"
            format = None
        d = dict(
            name=i,
            id=i,
            deletable=False,
            selectable=False,
            type=dtype,
            format=format
        )
        columns.append(d)
    data = df.to_dict('records')
    return data, columns, list(range(min(len(df), 3)))


@callback(
    Output('volcano-highlight-ids', 'data', allow_duplicate=True),
    Input('deseq-table', 'selected_rows'),
    State('volcano-highlight-ids', 'data'),
    State('dataset-dd', 'value'),
    State('comparison-dd', 'value'),
    prevent_initial_call='initial_duplicate'

)
def add_selected_rows(selected_rows, highlight_data, dataset_key, comp):
    df = get_deseq_result(dataset_key, comp)
    indices = df.iloc[selected_rows].index
    highlight_data["Selected"] = indices
    return highlight_data


@callback(
    Output("enrichment-graph", "figure"),
    Input("dataset-dd", "value"),
    Input("comparison-dd", "value"),
    Input("enrich-type-dd", "value"),
    Input("enrich-updown-dd", "value"),
    Input("mode-switch", "value"),
)
def create_enrich(dataset_key, comp, enrich_type, updown, switch):

    df = get_enrich_result(dataset_key, comp, enrich_type, updown)
    cu = UP_COLOR_LIGHT if switch else UP_COLOR_DARK
    cd = DOWN_COLOR_LIGHT if switch else DOWN_COLOR_DARK
    if df is not None:
        fig = enrichment_plot_from_cp_table(df, colorscale=[cu, cd])
    else:
        fig = empty_figure("No enrichment file found for dataset")
    if not switch:
        fig.update_layout(DARK_LAYOUT)
    else:
        fig.update_layout(LAYOUT)

    return fig


@callback(
    Output("correlation-graph", "figure"),
    Input("dataset-dd", "value"),
    Input("mode-switch", "value"),
)
def create_qc(dataset_key, switch):

    df = DASH_DATA[0][dataset_key]["qc"].get("correlation", None)
    cu = UP_COLOR_LIGHT if switch else UP_COLOR_DARK
    cd = DOWN_COLOR_LIGHT if switch else DOWN_COLOR_DARK
    tcol = "black" if switch else "white"
    bg_col = "white" if switch else None
    if df is not None:
        fig = pheatmap(df, colorscale=[cu, cd], tree_color=tcol, vertical_spacing=0, horizontal_spacing=0)
        if not switch:
            fig.update_layout(DARK_LAYOUT)
        else:
            fig.update_layout(LAYOUT)
            for p in ["", 4]:
                fig.add_shape(
                    type='rect',
                    xref=f'x{p}', yref=f'y{p}',
                    x0=-100, x1=100,
                    y0=-100, y1=100,  # row 1 (top)
                    fillcolor=bg_col,  # fully transparent
                    line=dict(width=0),
                    layer='below'
                )
        fig.update_yaxes(showline=False, showgrid=False, zeroline=False, col=2)
        fig.update_yaxes(showline=False, showgrid=False, zeroline=False, row=1, showticklabels=False)
        fig.update_xaxes(showline=False, showgrid=False, zeroline=False, col=2, showticklabels=False)
        fig.update_xaxes(showline=False, showgrid=False, zeroline=False, row=1)

    else:
        fig = empty_figure("No Correlation file found for dataset")
        if not switch:
            fig.update_layout(DARK_LAYOUT)
        else:
            fig.update_layout(LAYOUT)

    return fig


@callback(
    Output("pca-graph", "figure"),
    Input("dataset-dd", "value"),
    Input("mode-switch", "value"),
)
def create_qc(dataset_key, switch):
    df = DASH_DATA[0][dataset_key]["qc"].get("pca", None)
    cu = UP_COLOR_LIGHT if switch else UP_COLOR_DARK
    cd = DOWN_COLOR_LIGHT if switch else DOWN_COLOR_DARK
    if df is not None:
        colors = [cu, cd]
        colors = DEFAULT_PLOTLY_COLORS_LIST
        fig = sample_pca(df, colors=colors)
        fig.update_traces(marker=dict(size=10))
    else:
        fig = empty_figure("No PCA file found for dataset")
    if not switch:
        fig.update_layout(DARK_LAYOUT)
    else:
        fig.update_layout(LAYOUT)

    return fig





@callback(
    Output("comparison-dd", "value"),
    Output("comparison-dd", "options"),
    Input('dataset-dd', 'value'),
)
def update_comparisons(dataset_key):
    comps = list(DASH_DATA[0][dataset_key]["comparisons"].keys())
    selected = comps[0]
    return selected, comps



@callback(
    Output("enrich-updown-dd", "value"),
    Output("enrich-updown-dd", "options"),
    Input("comparison-dd", "value"),
    State('dataset-dd', 'value'),
)
def update_updown_dd(comp, dataset_key):
    condition = DASH_DATA[0][dataset_key]["comparisons"][comp]["condition"]
    baseline = DASH_DATA[0][dataset_key]["comparisons"][comp]["baseline"]
    return condition, [condition, baseline]



@callback(
    Output("add-name-dd", "value"),
    Output("add-name-dd", "options"),
    Input("comparison-dd", "value"),
    State('dataset-dd', 'value'),
    State('add-name-dd', 'value'),
)
def update_volcano_column_selections(comp, dataset_key, current_add_name):
    df = get_deseq_result(dataset_key, comp)
    columns = list(df.columns)
    sel = dash.no_update if current_add_name in columns else None
    return sel, columns


@callback(
    Output("volcano-graph", "figure"),
    Input('volcano-highlight-ids', 'data'),
    Input('add-name-dd', 'value'),
    Input('volcano-type-dd', 'value'),
    State('dataset-dd', 'value'),
    State('comparison-dd', 'value'),
    State('enrich-term', 'data'),
    State("mode-switch", "value"),

)
def create_volcano(highlight_data, name_col, volcano_type, dataset_key, comp, enrich_term, switch):
    highlight = {}
    for idx, (key, value) in enumerate(highlight_data.items()):
        highlight[key] = (DEFAULT_PLOTLY_COLORS[key], value)
    df = get_deseq_result(dataset_key, comp)
    config = DASH_DATA[0][dataset_key]["config"]
    if volcano_type == "Volcano":
        fig = volcano_from_deseq_result(
            df,
            highlight=highlight,
            name_col=name_col,
            lfc_cutoff=config["log2FCCutOff"],
            padj_cutoff=config["pAdjCutOff"],
            condition_name=DASH_DATA[0][dataset_key]["comparisons"][comp]["condition"],
            base_name=DASH_DATA[0][dataset_key]["comparisons"][comp]["baseline"],
            highlight_up_color=UP_COLOR_LIGHT if switch else UP_COLOR_DARK,
            highlight_down_color=DOWN_COLOR_LIGHT if switch else DOWN_COLOR_DARK,
            opacity=0.1 if switch else 0.25
        )
    else:
        fig = ma_from_deseq_result(
            df,
            name_col=name_col,
            highlight=highlight,
            lfc_cutoff=config["log2FCCutOff"],
            padj_cutoff=config["pAdjCutOff"],
            condition_name=DASH_DATA[0][dataset_key]["comparisons"][comp]["condition"],
            base_name=DASH_DATA[0][dataset_key]["comparisons"][comp]["baseline"],
            highlight_up_color=UP_COLOR_LIGHT if switch else UP_COLOR_DARK,
            highlight_down_color=DOWN_COLOR_LIGHT if switch else DOWN_COLOR_DARK,
        )
    if enrich_term:
        if volcano_type == "MA":
            sigcolor = DEFAULT_PLOTLY_COLORS["placeholder"]
            notsigcolor = "#00004a" if switch else "#f6f1e3"
            sig_add = {"marker": dict(color=sigcolor)}
            nsig_add = {"marker": dict(color=notsigcolor)}
        else:
            sig_add = {}
            nsig_add = {}
        fig.update_traces(selector=dict(name='not-enriched'), legendgroup=enrich_term, legendgrouptitle=dict(text=enrich_term), **nsig_add)
        fig.update_traces(selector=dict(name='enriched'), legendgroup=enrich_term, legendgrouptitle=dict(text=enrich_term), **sig_add)
    fig.update_traces(selector=dict(name='Selected'), legendgroup="Table", legendgrouptitle=dict(text="Table"))
    if not switch:
        fig.update_layout(DARK_LAYOUT)
    else:
        fig.update_layout(LAYOUT)

    return fig

@callback(
    Output('reset-enrich-selection', 'disabled', allow_duplicate=True),
    Output('volcano-highlight-ids', 'data', allow_duplicate=True),
    Input('reset-enrich-selection', 'n_clicks'),
    State('volcano-highlight-ids', 'data'),
    prevent_initial_call=True

)
def reset_enrich_selection(n_clicks, current_data):
    del current_data["enriched"]
    del current_data["not-enriched"]
    return True, current_data


@callback(
    Output('volcano-highlight-ids', 'data', allow_duplicate=True),
    Output('enrich-term', 'data', allow_duplicate=True),
    Output('reset-enrich-selection', 'disabled'),
    Input("enrichment-graph", "clickData"),
    State('volcano-highlight-ids', 'data'),
    State('dataset-dd', 'value'),
    State('comparison-dd', 'value'),
    State("enrich-type-dd", "value"),
    State("enrich-updown-dd", "value"),
    prevent_initial_call=True

)
def update_volcano_from_enrich(click_data, current_data, dataset_key, comp, enrich_type, updown):
    enrich = get_enrich_result(dataset_key, comp, enrich_type, updown)

    if click_data is not None:
        category = click_data["points"][0]["y"]
        sdf = enrich[enrich["Description"] == category].iloc[0]
        ids = set(sdf["geneID"].split("/"))
        if "universeGeneID" in sdf:
            ids2 = set(sdf["universeGeneID"].split("/"))
            ids2 = ids2 - ids2.intersection(ids)
            current_data["not-enriched"] = list(ids2)

        current_data["enriched"] = list(ids)
        print("CD", current_data)
        return current_data, category, False
    raise PreventUpdate




@callback(
    Output("volcano-graph", "figure", allow_duplicate=True),
    Input("mode-switch", "value"),
    State("volcano-graph", "figure"),
    State("gsea-graph", "figure"),
    State('volcano-type-dd', 'value'),
    prevent_initial_call='initial_duplicate'

)
def patch_all_figures_style(switch_value, f1, f2, volcano_type):
    upc = UP_COLOR_LIGHT if switch_value else UP_COLOR_DARK
    downc = DOWN_COLOR_LIGHT if switch_value else DOWN_COLOR_DARK
    font_color = "black" if switch_value else "white"
    patched_volcano = Patch()
    patched_volcano['layout']['font']['color'] = font_color
    patched_volcano['layout']['xaxis']['linecolor'] = font_color
    patched_volcano['layout']['yaxis']['linecolor'] = font_color

    patched_volcano["layout"]['plot_bgcolor']= LAYOUT["plot_bgcolor"] if switch_value else DARK_LAYOUT["plot_bgcolor"]
    if volcano_type == "Volcano":
        patched_volcano["layout"]["shapes"][0]["fillcolor"] = upc
        patched_volcano["layout"]["shapes"][0]["opacity"] = 0.1 if switch_value else 0.25
        patched_volcano["layout"]["annotations"][0]["font"]["color"] = upc
        patched_volcano["layout"]["shapes"][1]["fillcolor"] = downc
        patched_volcano["layout"]["shapes"][1]["opacity"] = 0.1 if switch_value else 0.25
        patched_volcano["layout"]["annotations"][1]["font"]["color"] = downc
    else:
        patched_volcano["data"][1]["marker"]["color"] = downc
        patched_volcano["data"][4]["marker"]["color"] = "#00004a" if switch_value else "#f6f1e3"

    return patched_volcano if f1 is not None else dash.no_update



