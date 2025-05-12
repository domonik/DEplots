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
import numpy as np



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
    comp = comps[0]
    sets = get_datasets_with_comp(comp)
    dataset_card = dbc.Col(
        dbc.Card(
            dbc.Row(
                [
                    dbc.Col(html.H5("Datasets"), width=2, className="d-flex"),
                    dbc.Col(
                        dcc.Dropdown(
                            comps, comp,
                            id="comparison-hash-dd",
                            clearable=False,
                            persistence=True,
                            persistence_type="session"

                        ),
                        width=3

                    ),
                    dbc.Col(
                        dcc.Dropdown(
                            id="dataset-compare-dd",
                            value=sets,
                            options=sets,
                            clearable=True,
                            multi=True,
                            persistence=True,
                            persistence_type="session",

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
                        dbc.Col([
                            html.H5("Multi DESeq table"),

                        ]
                            , width=3, align="center",
                        ),
                        dbc.Col(
                            html.Span(
                                [
                                    html.I(className="fas fa-xl fa-question-circle fa px-2",
                                           id="filter-tip"),
                                    dbc.Tooltip(
                                        "You can filter this table for up/down-regulated genes in the"
                                        " intersection of datasets. You can also click on the corresponding "
                                        "intersection in the upset "
                                        "plot below the table.",
                                        target="filter-tip"),
                                ],

                            ),
                            width=1,
                            className="d-flex align-items-center justify-content-end"

                        ),
                        dbc.Col(
                            [
                                dcc.Dropdown(
                                    id="filter-table-updown",
                                    clearable=True,
                                    persistence=True,
                                    persistence_type="session"

                                ),


                            ],
                            width=3

                        ),
                        dbc.Col(
                            [
                                dcc.Dropdown(
                                    options=["in", "exclusive in"],
                                    value="in",
                                    id="filter-table-operation",
                                    clearable=False,
                                    persistence=True,
                                    persistence_type="session"

                                ),

                            ],
                            width=1

                        ),
                        dbc.Col(
                            dcc.Dropdown(
                                id="filter-table-sets",

                                clearable=True,
                                multi=True,
                                persistence=True,
                                persistence_type="session",

                            ),
                            width=4
                        ),

                    ],
                    justify="between"
                ),

            ),
            dbc.Row(
                [

                    dbc.Col(
                        dash_table.DataTable(
                            id='datasets-table',
                            filter_action="custom",
                            filter_options={"case": "insensitive"},
                            sort_action="custom",
                            sort_mode="multi",
                            column_selectable="single",
                            row_selectable="multi",
                            row_deletable=False,
                            merge_duplicate_headers=True,
                            selected_columns=[],
                            selected_rows=[],
                            page_action="custom",

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
                            dbc.Col(html.H5(f"Upset {updown} regulated", id=f"{updown}-reg-header"), width=6),

                        ]

                    ),

                ),
                dbc.Col(dcc.Graph(id=f"upset-{updown}-graph",), width=12,

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
        padj_cutoff = DASH_DATA[0][datasets[0]]["config"]["pAdjCutOff"]
        if updown == "up":
            lfc_cutoff = DASH_DATA[0][datasets[0]]["config"]["log2FCCutOff"]
            barcolor = UP_COLOR_LIGHT if switch else UP_COLOR_DARK
            dot_colors = (dot_color, barcolor)

        else:
            lfc_cutoff = -DASH_DATA[0][datasets[0]]["config"]["log2FCCutOff"]
            barcolor = DOWN_COLOR_LIGHT if switch else DOWN_COLOR_DARK
            dot_colors = (dot_color, barcolor)
        fig = upset_plot_from_deseq(df, padj_cutoff=padj_cutoff, lfc_cutoff=lfc_cutoff, vertical_spacing=0, bar_color=barcolor, dot_colors=dot_colors, horizontal_spacing=0, mode=updown)
        fig.update_yaxes(fixedrange=True, row=2)
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
    Output("down-reg-header", "children"),
    Output("up-reg-header", "children"),
    Output("filter-table-updown", "options"),
    Output("filter-table-updown", "value"),
    Input("comparison-hash-dd", "value"),

)
def update_name_selection(comp):
    if comp is None:
        raise dash.exceptions.PreventUpdate
    datasets = get_datasets_with_comp(comp)
    options = [f"{col[0]} - {col[1]}" for col in DASH_DATA[1][comp].columns if col[0] not in datasets]
    data = DASH_DATA[0][datasets[0]]["comparisons"][comp]
    condition = data["condition"]
    baseline = data["baseline"]
    tableoptions = [condition, baseline]
    return options, options[0], f"Upset plot - {baseline}", f"Upset plot - {condition}", tableoptions, None


@callback(
    Output("filter-table-sets", "options"),
    Output("filter-table-sets", "value"),
    Input("dataset-compare-dd", "value"),
    State("filter-table-sets", "value"),

)
def update_set_selection(datasets, selected_sets):
    if datasets is None:
        return dash.no_update
    selected = None
    if selected_sets is not None:
        selected = selected_sets if all(ds in datasets for ds in selected_sets) else None
    return datasets, selected


@callback(
    Output("filter-table-sets", "value", allow_duplicate=True),
    Output("filter-table-updown", "value", allow_duplicate=True),
    Output("filter-table-operation", "value", allow_duplicate=True),
    Input("upset-up-graph", "clickData"),
    State("dataset-compare-dd", "value"),
    State("comparison-hash-dd", "value"),
    prevent_initial_call=True

)
def update_set_on_click(click_data, datasets, comp):
    if click_data is None:
        return dash.no_update
    s = click_data["points"][0]["x"]
    if isinstance(s, int):
        return dash.no_update
    else:
        sets = s.split("<br>")
    upreg = DASH_DATA[0][datasets[0]]["comparisons"][comp]["condition"]
    return sets, upreg, "exclusive in"

@callback(
    Output("filter-table-sets", "value", allow_duplicate=True),
    Output("filter-table-updown", "value", allow_duplicate=True),
    Output("filter-table-operation", "value", allow_duplicate=True),

    Input("upset-down-graph", "clickData"),
    State("dataset-compare-dd", "value"),
    State("comparison-hash-dd", "value"),
    prevent_initial_call=True

)
def update_set_on_click(click_data, datasets, comp):
    if click_data is None:
        return dash.no_update
    s = click_data["points"][0]["x"]
    if isinstance(s, int):
        return dash.no_update
    else:
        sets = s.split("<br>")
    upreg = DASH_DATA[0][datasets[0]]["comparisons"][comp]["baseline"]
    return sets, upreg, "exclusive in"

operators = [
    ['contains '],
    ['ge ', '>='],
    ['le ', '<='],
    ['lt ', '<'],
    ['gt ', '>'],
    ['ne ', '!='],
    ['eq ', '='],
    ['datestartswith ']
]


def split_filter_part(filter_part):
    for operator_type in operators:
        for operator in operator_type:
            if operator in filter_part:
                name_part, value_part = filter_part.split(operator, 1)
                print(name_part, value_part, "vname")
                name = name_part[name_part.find('{') + 1: name_part.rfind('}')]

                value_part = value_part.strip()
                v0 = value_part[0]
                if (v0 == value_part[-1] and v0 in ("'", '"', '`')):
                    value = value_part[1: -1].replace('\\' + v0, v0)
                else:
                    try:
                        value = float(value_part)
                    except ValueError:
                        value = value_part

                # word operators need spaces after them in the filter string,
                # but we don't want these later
                return name, operator_type[0].strip(), value

    return [None] * 3

@callback(
    Output("datasets-table", "columns"),
    Output("datasets-table", "data"),
    Output("datasets-table", "page_count"),
    Input("dataset-compare-dd", "value"),
    Input("filter-table-sets", "value"),
    Input("filter-table-updown", "value"),
    Input("filter-table-operation", "value"),
    Input('datasets-table', "page_current"),
    Input('datasets-table', "page_size"),
    Input('datasets-table', 'sort_by'),
    Input('datasets-table', 'filter_query'),
    State("comparison-hash-dd", "value"),
)
def update_datasets_table(datasets, filter_set, filter_ud, filter_op, current_page, page_size, sort_by, filter_query, comp):

    if datasets is None:
        return None, None
    idx = pd.IndexSlice
    df = DASH_DATA[1][comp]
    df["id"] = range(len(df))

    df = df.loc[:, idx[["id", "Name", "Additional Data"] + datasets, :]]
    exclusive = filter_op == "exclusive in"
    if filter_set is not None and len(filter_set) == 0:
        filter_set = None
    if filter_ud and filter_set:
        upreg = DASH_DATA[0][datasets[0]]["comparisons"][comp]["condition"] == filter_ud
        padj_cutoff = DASH_DATA[0][datasets[0]]["config"]["pAdjCutOff"]
        lfc_cutoff = DASH_DATA[0][datasets[0]]["config"]["log2FCCutOff"]
        idx = pd.IndexSlice
        others = list(set(datasets) - set(filter_set))
        lfc_idx = idx[datasets, ('log2FoldChange')]
        pad_idx = idx[datasets, ('padj')]
        lfc_df = df.loc[:, lfc_idx].droplevel(1, axis=1)
        padj_df = df.loc[:, pad_idx].droplevel(1, axis=1)
        if upreg:
            set_idx = np.all((lfc_df.loc[:, filter_set] >= lfc_cutoff) & (padj_df.loc[:, filter_set] <= padj_cutoff), axis=1)
            if exclusive:
                other_idx = np.any((lfc_df.loc[:, others] >= lfc_cutoff) & (padj_df.loc[:, others] <= padj_cutoff), axis=1)
                set_idx = (set_idx & ~other_idx)

        else:
            lfc_cutoff = -lfc_cutoff
            set_idx = np.all((lfc_df.loc[:, filter_set] <= lfc_cutoff) & (padj_df.loc[:, filter_set] <= padj_cutoff),
                             axis=1)
            if exclusive:
                other_idx = np.any((lfc_df.loc[:, others] <= lfc_cutoff) & (padj_df.loc[:, others] <= padj_cutoff), axis=1)
                set_idx = (set_idx & ~other_idx)
        df = df[set_idx]
    if filter_query:
        filtering_expressions = filter_query.split(' && ')
        for filter_part in filtering_expressions:
            col_name, operator, filter_value = split_filter_part(filter_part)
            col_name = tuple(col_name.split("_"))


            if operator in ('eq', 'ne', 'lt', 'le', 'gt', 'ge'):
                # these operators match pandas series operator method names
                df = df.loc[getattr(df[col_name], operator)(filter_value)]
            elif operator == 'contains':
                df = df.loc[df[col_name].astype(str).str.contains(filter_value, case=False) == True]
            elif operator == 'datestartswith':
                # this is a simplification of the front-end filtering logic,
                # only works with complete fields in standard format
                df = df.loc[df[col_name].str.startswith(filter_value)]
    if sort_by:
        sort_by_filtered = []
        for i, col in enumerate(sort_by):
            t = tuple(col["column_id"].split("_"))
            if t in df.columns:
                sort_by_filtered.append( sort_by[i] |{"column_id": t})
        sort_by = sort_by_filtered
        if len(sort_by):
            sort_vals = [col['column_id'] for col in sort_by if col["column_id"]]
            print("sort vals", sort_vals)
            df = df.sort_values(
                sort_vals,
                ascending=[
                    col['direction'] == 'asc'
                    for col in sort_by
                ],
                inplace=False
            )
    page_count = (df.shape[0] // page_size) + 1
    df = df.iloc[current_page*page_size:(current_page+ 1)*page_size]
    columns = [
        {
            "name": i,
            "id": "_".join(i) if i[0] != "id" else "id",
            "deletable": False,
            "selectable": False,
            "format": Format(precision=4),
            "type": "numeric" if is_numeric_dtype(df[i]) else "text"
        } for i in df.columns if i[0] != "id"
    ]
    data = [{"_".join(col) if col[0] != "id" else "id": val for col, val in row.items()} for row in df.to_dict('records')]
    return columns, data, page_count


@callback(
    Output("dataset-compare-dd", "options"),
    Output("dataset-compare-dd", "value"),
    Input("comparison-hash-dd", "value"),
    State("dataset-compare-dd", "value"),

)
def update_selectable_datasets(comp, dd_state):
    sets = get_datasets_with_comp(comp)
    if dd_state is not None:
        selected = dd_state if all([v in sets for v in dd_state]) else sets
    else:
        selected = sets
    return sets, selected



@callback(
    Output("gene-line-graph", "figure"),
    Input("datasets-table", "derived_virtual_selected_row_ids"),
    Input("dataset-compare-dd", "value"),
    Input("mode-switch", "value"),
    Input("plot-hover-name-dd", "value"),

    State("comparison-hash-dd", "value"),

)
def update_line_plot(sel_rows, datasets, switch, legend_name, comp):
    plot = False
    print(sel_rows, "sel_rows", dash.ctx.triggered_prop_ids)
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