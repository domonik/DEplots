from DEplots.volcano import volcano_from_deseq_result
from DEplots.enrichment import enrichment_plot_from_cp_table, empty_figure, plot_gsea
import dash
from dash import html, clientside_callback, Input, Output, Dash, dcc, dash_table, State, Patch
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import pandas as pd
import os
from DEplots.dashboard import get_data
from dash.dash_table.Format import Format, Scheme, Trim
from pandas.api.types import is_numeric_dtype
import plotly.express as px


FILEDIR = os.path.dirname(os.path.abspath(__file__))
ASSETS_DIR = os.path.join(FILEDIR, "assets")
assert os.path.exists(ASSETS_DIR)

DEFAULT_PLOTLY_COLORS = {
    "not-enriched": "#344A9A",
    "enriched": "#00a082",
    "Selected": "#8f6b30",
    "placeholder": "#ffe863",
    "placeholder2": "#f5c2ed"

}
DEFAULT_PLOTLY_COLORS_LIST = list(DEFAULT_PLOTLY_COLORS.values()) + px.colors.DEFAULT_PLOTLY_COLORS



LAYOUT = {
    "template": "plotly_white",
    'paper_bgcolor': 'rgba(0,0,0,0)',
    'plot_bgcolor': 'rgba(0,0,0,0)',
    "font": {"color": "black", "size": 16},
    "xaxis": {"showline": True, "mirror": True, "linecolor": "black"},
    "yaxis": {"showline": True, "mirror": True, "linecolor": "black"},
    "margin": {"b": 10, "t": 10}
}

DARK_LAYOUT = {
    "template": "plotly_white",
    'paper_bgcolor': 'rgba(0,0,0,0)',
    'plot_bgcolor': 'rgba(0,0,0,0)',
    "font": {"color": "white", "size": 16},
    "xaxis": {"showline": True, "mirror": True, "linecolor": "white"},
    "yaxis": {"showline": True, "mirror": True, "linecolor": "white"},
    "margin": {"b": 10, "t": 10}
}


UP_COLOR_LIGHT = "#00a082"
UP_COLOR_DARK = "#00a082"
DOWN_COLOR_LIGHT = "#344A9A"
DOWN_COLOR_DARK = "#f5c2ed"

dbc_css = "https://cdn.jsdelivr.net/gh/AnnMarieW/dash-bootstrap-templates/dbc.min.css"

app = Dash(
    __name__,
    external_stylesheets=["custom.css", dbc.icons.FONT_AWESOME, dbc_css],
    assets_folder=ASSETS_DIR,
    prevent_initial_callbacks=True
)


def get_deseq_result(dataset_key, comp):
    df = DASH_DATA[dataset_key]["comparisons"][comp]["deseq"]
    return df


def get_enrich_result(dataset_key, comp, enrich: str = "GO", updown: str = "up-regulated"):
    df = DASH_DATA[dataset_key]["comparisons"][comp]["enrich"][enrich][updown]
    return df

def get_gsea_result(dataset_key, comp):
    df = DASH_DATA[dataset_key]["comparisons"][comp]["gsea"]["df"]
    plot_data = DASH_DATA[dataset_key]["comparisons"][comp]["gsea"]["plot_data"]
    return df, plot_data



color_mode_switch =  dbc.Row(
    [
        dbc.Col(html.I(className="fa fa-xl fa-moon", ), className="d-flex justify-content-end align-items-center", width=3, md=1),
        dbc.Col(dbc.Switch(id="mode-switch", className="ms-3 fs-4", persistence=True), className="d-flex justify-content-center align-items-center", width=3, md=1),
        dbc.Col(html.I(className="fa fa-xl fa-sun"), className="d-flex justify-content-start align-items-center", width=3, md=1 ),
    ], justify="end", className="w-100"
)
PLOTLY_LOGO = "https://images.plot.ly/logo/new-branding/plotly-logomark.png"


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

navbar = dbc.Navbar(
    dbc.Container(
        [
            dbc.Row(
                [
                    dbc.Col(
                        dbc.Row(
                            [
                                dbc.Col(html.Img(src="https://cd.uni-freiburg.de/wp-content/uploads/2022/09/ufr-logo-white-2.png", height="50px")),
                            ],
                            align="center",
                            className="g-0",
                        ),
                        width=3, className="d-md-flex d-none"
                    ),
                    dbc.Col(html.H2("DESeq Explorer", className="text-md-center"), width=6, align="center",),

                    dbc.Col(
                        color_mode_switch,
                        width=3, className="d-flex align-items-center"
                    )
                ],
                className="w-100 ", justify="between"

            ),



        ],
        fluid=True,
        className="dbc text-light"

    ),
    color="var(--bs-ufr-navbar)",
    className="ufr-navbar shadow w-100", style={"position": "fixed", "z-index": "10"}
)


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


def get_layout(dash_data):
    layout = html.Div([
        navbar,
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
                                                dbc.Col(html.H5("DEseq Volcano"), width=6),
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
                            width=12, md=6
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


                                                    ), width=3, className="d-flex align-items-center ",
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
                            width=12, md=6
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

            ],
            fluid=True,
            className="dbc",
            style={"padding-top": "65px"}

        ),
        html.Footer(
            dbc.Container(
                [
                    dbc.Row(
                        [
                            dbc.Col(width=3),
                            dbc.Col(width=3),
                            dbc.Col(
                                html.Ul(
                                    [
                                        html.Li(html.A(html.I(className="fa-brands fa-2xl fa-github"), target="_blank", href="https://github.com/domonik/DEplots", className="text-light")),
                                        html.Li(html.A(html.I(className="fa-solid fa-2xl fa-envelope"), target="_blank", href="mailto:rabsch@informatik.uni-freiburg.de", className="text-light")),
                                        html.Li(html.A(html.I(className="fa-brands fa-2xl fa-linkedin"), target="_blank", href="https://www.linkedin.com/in/dominik-rabsch/", className="text-light"))
                                    ],
                                    className="icon-list d-flex align-items-end justify-content-end text-light"
                                ),
                                width=3,
                                align="end"
                            ),
                        ],
                        className="w-100 py-4", justify="between"
                    )
                ],
                fluid=True
            ),
            className="ufr-navbar text-light mt-4", style={"z-index": "20"}
        )
    ])
    return layout


@app.callback(
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


@app.callback(
    Output('gsea-graph', 'figure'),
    Input('gsea-table', 'selected_rows'),
    Input('mode-switch', 'value'),

    State('dataset-dd', 'value'),
    State('comparison-dd', 'value'),
    prevent_initial_call=False
)
def get_gsea_plot(selected_rows, switch, dataset_key, comp):
    df, plot_data = get_gsea_result(dataset_key, comp)
    descs = df["Description"].iloc[selected_rows]
    cond = DASH_DATA[dataset_key]["comparisons"][comp]["condition"]
    baseline = DASH_DATA[dataset_key]["comparisons"][comp]["baseline"]
    if df is not None and len(df) > 0:
        fig = plot_gsea(plot_data, descs=descs, colors=DEFAULT_PLOTLY_COLORS_LIST, show_zero_lfc=True, condition_name=cond, base_name=baseline, gene_list_name="log2FC")
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




@app.callback(
    Output("deseq-table", "data"),
    Output("deseq-table", "columns"),
    Output('volcano-highlight-ids', 'data'),
    Input('dataset-dd', 'value'),
    Input('comparison-dd', 'value'),
    prevent_initial_call=False
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


@app.callback(
    Output("gsea-table", "data"),
    Output("gsea-table", "columns"),
    Output('gsea-table', 'selected_rows'),
    Input('dataset-dd', 'value'),
    Input('comparison-dd', 'value'),
    prevent_initial_call=False
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


@app.callback(
    Output('volcano-highlight-ids', 'data', allow_duplicate=True),
    Input('deseq-table', 'selected_rows'),
    State('volcano-highlight-ids', 'data'),
    State('dataset-dd', 'value'),
    State('comparison-dd', 'value'),

)
def add_selected_rows(selected_rows, highlight_data, dataset_key, comp):
    df = get_deseq_result(dataset_key, comp)
    indices = df.iloc[selected_rows].index
    highlight_data["Selected"] = indices
    return highlight_data


@app.callback(
    Output("enrichment-graph", "figure"),
    Input("dataset-dd", "value"),
    Input("comparison-dd", "value"),
    Input("enrich-type-dd", "value"),
    Input("enrich-updown-dd", "value"),
    Input("mode-switch", "value"),

    prevent_initial_call=False
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


@app.callback(
    Output("comparison-dd", "value"),
    Output("comparison-dd", "options"),
    Input('dataset-dd', 'value'),
    prevent_initial_call=False

)
def update_comparisons(dataset_key):
    comps = list(DASH_DATA[dataset_key]["comparisons"].keys())
    selected = comps[0]
    return selected, comps

@app.callback(
    Output("enrich-updown-dd", "value"),
    Output("enrich-updown-dd", "options"),
    Input("comparison-dd", "value"),
    State('dataset-dd', 'value'),
)
def update_updown_dd(comp, dataset_key):
    condition = DASH_DATA[dataset_key]["comparisons"][comp]["condition"]
    baseline = DASH_DATA[dataset_key]["comparisons"][comp]["baseline"]
    return condition, [condition, baseline]

@app.callback(
    Output("add-name-dd", "value"),
    Output("add-name-dd", "options"),
    Input("comparison-dd", "value"),
    State('dataset-dd', 'value'),
    State('add-name-dd', 'value'),
)
def update_volcano_column_selections(comp, dataset_key, current_add_name):
    df = DASH_DATA[dataset_key]["comparisons"][comp]["deseq"]
    columns = list(df.columns)
    sel = dash.no_update if current_add_name in columns else None
    return sel, columns

@app.callback(
    Output("volcano-graph", "figure"),
    Input('volcano-highlight-ids', 'data'),
    Input('add-name-dd', 'value'),
    State('dataset-dd', 'value'),
    State('comparison-dd', 'value'),
    State('enrich-term', 'data'),
    State("mode-switch", "value"),

)
def create_volcano(highlight_data, name_col, dataset_key, comp, enrich_term, switch):
    highlight = {}
    for idx, (key, value) in enumerate(highlight_data.items()):
        highlight[key] = (DEFAULT_PLOTLY_COLORS[key], value)
    df = get_deseq_result(dataset_key, comp)
    config = DASH_DATA[dataset_key]["config"]

    fig = volcano_from_deseq_result(
        df,
        highlight=highlight,
        name_col=name_col,
        lfc_cut_off=config["log2FCCutOff"],
        padj_cutoff=config["pAdjCutOff"],
        condition_name=DASH_DATA[dataset_key]["comparisons"][comp]["condition"],
        base_name=DASH_DATA[dataset_key]["comparisons"][comp]["baseline"],
        highlight_up_color=UP_COLOR_LIGHT if switch else UP_COLOR_DARK,
        highlight_down_color=DOWN_COLOR_LIGHT if switch else DOWN_COLOR_DARK,
        opacity=0.1 if switch else 0.25
    )
    if enrich_term:
        print(enrich_term)
        fig.update_traces(selector=dict(name='not-enriched'), legendgroup=enrich_term, legendgrouptitle=dict(text=enrich_term))
        fig.update_traces(selector=dict(name='enriched'), legendgroup=enrich_term, legendgrouptitle=dict(text=enrich_term))
    fig.update_traces(selector=dict(name='Selected'), legendgroup="Table", legendgrouptitle=dict(text="Table"))
    if not switch:
        fig.update_layout(DARK_LAYOUT)
    else:
        fig.update_layout(LAYOUT)

    return fig


@app.callback(
    Output('reset-enrich-selection', 'disabled', allow_duplicate=True),
    Output('volcano-highlight-ids', 'data', allow_duplicate=True),
    Input('reset-enrich-selection', 'n_clicks'),
    State('volcano-highlight-ids', 'data'),

)
def reset_enrich_selection(n_clicks, current_data):
    del current_data["enriched"]
    del current_data["not-enriched"]
    return True, current_data


@app.callback(
    Output('volcano-highlight-ids', 'data', allow_duplicate=True),
    Output('enrich-term', 'data', allow_duplicate=True),
    Output('reset-enrich-selection', 'disabled'),
    Input("enrichment-graph", "clickData"),
    State('volcano-highlight-ids', 'data'),
    State('dataset-dd', 'value'),
    State('comparison-dd', 'value'),
    State("enrich-type-dd", "value"),
    State("enrich-updown-dd", "value"),

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




@app.callback(
    Output("volcano-graph", "figure", allow_duplicate=True),
    Input("mode-switch", "value"),
    State("volcano-graph", "figure"),
    State("gsea-graph", "figure"),

)
def patch_all_figures_style(switch_value, f1, f2):
    upc = UP_COLOR_LIGHT if switch_value else UP_COLOR_DARK
    downc = DOWN_COLOR_LIGHT if switch_value else DOWN_COLOR_DARK
    font_color = "black" if switch_value else "white"
    patched_volcano = Patch()
    patched_volcano['layout']['font']['color'] = font_color
    patched_volcano['layout']['xaxis']['linecolor'] = font_color
    patched_volcano['layout']['yaxis']['linecolor'] = font_color
    patched_volcano["layout"]["shapes"][0]["fillcolor"] = upc
    patched_volcano["layout"]["shapes"][0]["opacity"] = 0.1 if switch_value else 0.25
    patched_volcano["layout"]["annotations"][0]["font"]["color"] = upc

    patched_volcano["layout"]["shapes"][1]["fillcolor"] = downc
    patched_volcano["layout"]["shapes"][1]["opacity"] = 0.1 if switch_value else 0.25
    patched_volcano["layout"]["annotations"][1]["font"]["color"] = downc

    return patched_volcano if f1 is not None else dash.no_update



clientside_callback(
    """
    (switchOn) => {
       document.documentElement.setAttribute("data-bs-theme", switchOn ? "light" : "dark"); 
       return window.dash_clientside.no_update
    }
    """,
    Output("mode-switch", "id"),
    Input("mode-switch", "value"),
)


def cli_wrapper(
        config_file: str = None,
        run_dir: str = None,
        debug: bool = False,
        port: int = 8080,
        host: str = "127.0.0.1",
        processes: int = 1
):
    global DASH_DATA
    DASH_DATA = get_data(config_file, run_dir)

    app.layout = get_layout(dash_data=DASH_DATA)
    app.run(debug=debug, port=port, host=host, processes=processes, threaded=False)


def _cli_wrapper(args):
    cli_wrapper(args.config, args.run_dir, args.debug, args.port, args.host, args.processes)



if __name__ == '__main__':
    config_file = "/home/rabsch/PythonProjects/DEPlots/testData/config.yaml"
    rd = "/home/rabsch/PythonProjects/RlocSeq/Pipeline/RUNS/"

    cli_wrapper(config_file=config_file, run_dir=rd, debug=True, processes=3)