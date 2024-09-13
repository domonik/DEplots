import dash
import gffutils
import pandas as pd
from dash import Input, Output, State, html, dcc, callback, dash_table
import dash_bootstrap_components as dbc
from dash.dash_table.Format import Format, Scheme
from pandas.api.types import is_numeric_dtype
from DEplots.runComparison import plot_gene_among_conditions, upset_plot_from_deseq
from DEplots.dashboard import DEFAULT_PLOTLY_COLORS, DEFAULT_PLOTLY_COLORS_LIST, LAYOUT, DARK_LAYOUT, UP_COLOR_LIGHT, \
    UP_COLOR_DARK, DOWN_COLOR_LIGHT, DOWN_COLOR_DARK, COVERAGE_DATA, COVERAGE_DESIGN, GFF, LINE_MAPPING
import numpy as np
from DEplots.readCount import plot_precomputed_coverage, GFF3_TYPE_TO_COLOR, css_color_to_rgba
import time
from DEplots.gff_helper import GFF_COLNAMES
import math

FEATURE_COLORS = GFF3_TYPE_TO_COLOR | {
    "gene": "rgb(219, 219, 219)",
    "mRNA": "rgb(219, 219, 219)",
    "exon": DOWN_COLOR_DARK,
    "CDS": UP_COLOR_DARK,
    "start_codon": UP_COLOR_LIGHT,
    "stop_codon": UP_COLOR_LIGHT,
    "five_prime_UTR": DOWN_COLOR_LIGHT,
    "5UTR": DOWN_COLOR_LIGHT,
    "5'UTR": DOWN_COLOR_LIGHT,
    "three_prime_UTR": DOWN_COLOR_LIGHT,
    "3UTR": DOWN_COLOR_LIGHT,
    "3'UTR": DOWN_COLOR_LIGHT,
    "ncRNA": UP_COLOR_DARK,
    # "rRNA": px.colors.qualitative.Light24[9],
    # "tRNA": px.colors.qualitative.Light24[10],
    # "repeat_region": px.colors.qualitative.Light24[11],
    "default": "#7CB9E8"
}
if COVERAGE_DATA is not None:
    dash.register_page(__name__, path='/coverage', name="Read Coverage")
    TRACE_COLORS = {trace: DEFAULT_PLOTLY_COLORS_LIST[idx] for idx, trace in
                    enumerate(COVERAGE_DESIGN.Treatment.unique())}
SETTINGS_ROW = "px-2"
SETTINGS_COL = "py-2 d-flex align-items-center justify-content-between"
SETTINGS_LABEL = {"className": "me-3", "style": {"min-width": "30%"}}


def coverage_settings_card():
    contigs = list(COVERAGE_DATA.keys())
    coverage_plot_card = dbc.Col(
        dbc.Card(
            [
                # dbc.CardHeader(
                #     dbc.Row(
                #         [
                #             dbc.Col(html.H5("Settings"), width=6),
                #
                #         ]
                #
                #     ),
                #
                # ),
                dbc.Col(
                    [
                        dbc.Row(
                            [
                                dbc.Col([
                                    html.Span("Sequence", **SETTINGS_LABEL, ),
                                    dcc.Dropdown(id="contig", options=contigs, value=contigs[0], clearable=False,
                                                 style={"width": "100%"})
                                ], width=12, sm=6, lg=3, className=SETTINGS_COL),
                                dbc.Col([
                                    html.Span("Start", **SETTINGS_LABEL),
                                    dbc.Input(id="coverage-start", type="number", value=0, min=0, step=1, ),

                                ], width=12, sm=6, lg=3, className=SETTINGS_COL),
                                dbc.Col([
                                    html.Span("End", **SETTINGS_LABEL),
                                    dbc.Input(id="coverage-end", type="number", value=2000, min=0, step=1, ),
                                ], width=12, sm=6, lg=3, className=SETTINGS_COL),
                                dbc.Col([
                                    html.Span("Autorange Y", **SETTINGS_LABEL),
                                    dbc.Switch(id="autorange-y", value=True, className="fs-5")
                                ], width=12, sm=6, lg=3, className=SETTINGS_COL),

                            ],
                            className=SETTINGS_ROW

                        ),

                        dbc.Row(
                            [

                            ],
                            className=SETTINGS_ROW

                        ),

                    ],
                    className="p-2"

                ),

            ],
            className="shadow",
        ),
        width=12, md=12
    )

    return coverage_plot_card


def coverage_card():
    coverage_plot_card = dbc.Col(
        dbc.Card(
            [
                dbc.CardHeader(
                    dbc.Row(
                        [
                            dbc.Col(html.H5("Coverage"), width=4),
                            dbc.Col(html.Span("Trace"), width=1,
                                    className="d-flex align-items-center justify-content-end"),
                            dbc.Col(dcc.Dropdown(
                                options=list(COVERAGE_DESIGN.Treatment.unique()),
                                value=list(COVERAGE_DESIGN.Treatment.unique())[0],
                                style={"width": "100%"},
                                id="trace-color-dd",
                                clearable=False

                            ), width=3, className="d-flex align-items-center"),
                            dbc.Col(
                                dbc.Input(
                                    id="trace-color",
                                    type="color",
                                    style={
                                        "width": "10rem",
                                        "height": "2rem",

                                    }
                                ),
                                width=1,
                                className="d-flex align-items-center"

                            )

                        ]

                    ),

                ),
                dbc.Col(
                    dcc.Graph(
                        id="coverage-graph",
                        style={"resize": "vertical", "overflow": "auto"},
                        config={
                            "scrollZoom": True,
                            'modeBarButtonsToRemove': ['zoom', 'autoscale', "resetScale2d"]

                        }
                    ), width=12,

                ),

            ],
            className="shadow",
        ),
        width=12
    )

    return coverage_plot_card


TYPE_MAP = {
               ftype: "text" for ftype in GFF_COLNAMES if ftype not in ["start", "end"]
           } | {"start": "numeric", "end": "numeric"}


def gff_container():
    if isinstance(GFF, str):
        data = gffutils.FeatureDB(GFF)
        cs = data.all_features()
        cols = GFF_COLNAMES

        data = [{colname: getattr(row, colname) for colname in GFF_COLNAMES if colname != "attributes"} | {
            "attributes": str(row).split("\t")[-1], "id": row.id} for row in cs]
        cols = [{"name": "id", "id": "id"}] + [{"name": i, "id": i, "type": TYPE_MAP[i]} for i in cols if i != "id"]
        print(cols)

    elif isinstance(GFF, pd.DataFrame):
        data = GFF
        data["id"] = data.index
        data = data.to_dict('records')
        cols = [{"name": i, "id": i} for i in GFF.columns if i != "id"]
    else:
        raise NotImplementedError()
    gff = dbc.Col(
        dbc.Card(
            [
                dbc.CardHeader(
                    dbc.Row(
                        [
                            dbc.Col(html.H5("GFF"), width=6),
                            dbc.Col(
                                [
                                    html.Div(html.Div(className="gff-legend2"), className="gff-legend"),
                                    html.Span("Feature visible", className="ms-3"),

                                ], width=6, className="d-flex align-items-center justify-content-end"),
                        ]

                    ),

                ),
                dbc.Col(
                    dash_table.DataTable(
                        id='gff-table',
                        columns=cols,
                        data=data,
                        editable=True,
                        filter_action="native",
                        filter_options={"case": "insensitive"},
                        sort_action="native",
                        sort_mode="multi",
                        column_selectable=False,
                        merge_duplicate_headers=True,
                        row_selectable=False,
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
                            # {
                            #     'if': {
                            #         'column_type': 'text'  # 'text' | 'any' | 'datetime' | 'numeric'
                            #     },
                            #     'textAlign': 'left'
                            # },

                        ],
                        style_cell_conditional=[
                            {
                                'if': {'column_id': 'attributes'},
                                'textAlign': 'left',

                            }
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
                    width=12,
                    style={"overflow": "auto", 'backgroundColor': 'var(--bs-primary-bg)'},
                    className="p-2 gff-table"

                ),

            ],
            className="shadow",
        ),
        width=12, lg=6,
    )
    return gff


def visible_gff_container():
    if isinstance(GFF, str):
        cols = [{"name": "id", "id": "id"}] + [{"name": i, "id": i} for i in GFF_COLNAMES if i != "id"]
    elif isinstance(GFF, pd.DataFrame):
        cols = [{"name": i, "id": i} for i in GFF.columns if i != "id"]
    else:
        raise NotImplementedError
    gff = dbc.Col(
        dbc.Card(
            [
                dbc.CardHeader(
                    dbc.Row(
                        [
                            dbc.Col(html.H5("Visible"), width=6),
                            dbc.Col(
                                [
                                    html.Div(html.Div(className="gff-legend2"), className="gff-legend"),
                                    html.Span("Feature visible", className="ms-3"),

                                ], width=6, className="d-flex align-items-center justify-content-end"),
                        ]

                    ),

                ),
                dbc.Col(
                    dash_table.DataTable(
                        id='visible-gff-table',
                        columns=cols,
                        data=None,
                        editable=True,
                        sort_action="native",
                        sort_mode="multi",
                        column_selectable=False,
                        merge_duplicate_headers=True,
                        row_selectable=False,
                        row_deletable=False,
                        selected_columns=[],
                        selected_rows=[],
                        page_action="native",
                        page_current=0,
                        page_size=10,
                        style_as_list_view=True,
                        style_data_conditional=[
                            {
                                'if': {'row_index': 'even'},
                                'backgroundColor': "rgba(var(--bs-primary-rgb), 0.2)",
                            },
                            {
                                'if': {'row_index': 'odd'},
                                'backgroundColor': "rgba(var(--bs-primary-rgb), 0.3)",
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
                        style_data={'border': 'none !important'},
                        style_cell_conditional=[
                            {
                                'if': {'column_id': 'attributes'},
                                'textAlign': 'left',

                            }
                        ],

                    ),
                    width=12,
                    style={"overflow": "auto", 'backgroundColor': 'var(--bs-primary-bg)'},
                    className="p-2"

                ),

            ],
            className="shadow",
        ),
        width=12, lg=6,
    )
    return gff


def get_layout():
    lout = html.Div(
        [
            dcc.Store(id="x-axis-range", data={"xaxis.range[0]": 1, "xaxis.range[1]": 2000, "yaxis2.range[0]": -0.5,
                                               "y2axis.range[1]": 5.5}),
            dcc.Store(id="internal-axis-range", data=[0, 2000]),
            dcc.Store(id="placeholder1", data=None),
            dcc.Store(id="current-contig", data=None),
            dcc.Store(id="trace-colors", data=TRACE_COLORS),
            dbc.Container(
                [
                    dbc.Row(
                        coverage_card(),
                        className="py-1"

                    ),
                    dbc.Row(
                        coverage_settings_card(),
                    ),
                    dbc.Row(
                        [

                            gff_container(),
                            visible_gff_container(),
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


def _is_zoom_relayout(old_rlout_data, new_rlout_data):
    old_winsize = old_rlout_data["xaxis.range[1]"] - old_rlout_data["xaxis.range[0]"]
    new_winsize = new_rlout_data["xaxis.range[1]"] - new_rlout_data["xaxis.range[0]"]
    if np.isclose(old_winsize, new_winsize):  # only shifts window
        return False
    else:
        return True


@callback(
    Output("trace-color", "value"),
    Input("trace-color-dd", "value"),
    State("trace-colors", "data"),
)
def update_box_color(trace, colors):
    return colors[trace]


@callback(
    Output("x-axis-range", "data", allow_duplicate=True),
    Output("coverage-start", "value", allow_duplicate=True),
    Output("coverage-end", "value", allow_duplicate=True),
    Input('coverage-graph', 'relayoutData'),
    State("x-axis-range", 'data'),
    State("internal-axis-range", "data"),

    prevent_initial_call=True
)
def detect_relayout(new_lout, old_lout, internal_xaxis):
    if new_lout is None:
        raise dash.exceptions.PreventUpdate
    if "xaxis.range[0]" in new_lout:
        start = max(int(new_lout["xaxis.range[0]"]), 0)
        end = int(new_lout["xaxis.range[1]"])
        winsize = end - start
        add_win = winsize * 0.25

        is_zoom = _is_zoom_relayout(old_lout, new_lout)
        if is_zoom:
            print("zoomed")
            print(old_lout)
            print(new_lout)
            new_lout["xaxis.range[0]"] = -1  # This forces the figure to replot

            return new_lout, start, end
        else:  # now we need to determine if we moved outside the window_range
            if new_lout["xaxis.range[1]"] > internal_xaxis[1] - add_win or new_lout["xaxis.range[0]"] < internal_xaxis[
                0] + add_win:
                print("scrolled out of window")
                new_lout["xaxis.range[0]"] = -1  # This forces the figure to replot
                print(new_lout)
                return new_lout, start, end
        return new_lout, start, end

    raise dash.exceptions.PreventUpdate


def _calc_internal_params(display_start, display_end):
    winsize = display_end - display_start
    add_win = winsize * 0.5
    start = display_start - add_win
    end = display_end + add_win
    step = max(1, int((end - start) * 0.005))
    return [start, end], step


@callback(
    Output("coverage-start", "value", allow_duplicate=True),
    Output("coverage-end", "value", allow_duplicate=True),
    Output("contig", "value", allow_duplicate=True),
    Input("gff-table", "active_cell"),
    prevent_initial_call=True
)
def click_data_table(active_cell):
    if active_cell is None:
        raise dash.exceptions.PreventUpdate
    if isinstance(GFF, pd.DataFrame):

        data = GFF.iloc[active_cell["row_id"]]
    elif isinstance(GFF, str):
        db = gffutils.FeatureDB(GFF)
        data = db[active_cell["row_id"]]
    else:
        raise NotImplementedError

    start = max(data.start - 250, 0)
    contig = data.seqid
    end = min(COVERAGE_DATA[contig]["+"].shape[-1] - 1, data.end + 250)

    return start, end, contig


@callback(
    Output("coverage-start", "value"),
    Output("coverage-end", "value"),
    Output("coverage-end", "max"),
    Input("contig", "value"),
    State("coverage-start", "value"),
    State("coverage-end", "value"),
)
def update_via_contig(contig, start, end):
    size = COVERAGE_DATA[contig]["+"].shape[-1]
    print(size, start, end)
    if end is None or end >= size:
        start = 0
        end = min(COVERAGE_DATA[contig]["+"].shape[-1], 2000)

    return start, end, size


def check_if_update_necessary(start, end, display_start, display_end, old_contig, contig):
    if dash.ctx.triggered_id != "mode-switch" and old_contig == contig:
        if any(i is None for i in (start, end)):
            raise dash.exceptions.PreventUpdate
        if int(display_start) == int(start) and int(display_end) == int(end):
            raise dash.exceptions.PreventUpdate
        if start >= end:
            raise dash.exceptions.PreventUpdate


@callback(
    Output("visible-gff-table", "data"),
    Input("coverage-start", "value"),
    Input("coverage-end", "value"),
    State("contig", "value"),
    State("current-contig", "data"),
    State("x-axis-range", "data"),

)
def update_visible_table(start, end, contig, old_contig, axis_range):
    display_start = axis_range["xaxis.range[0]"]
    display_end = axis_range["xaxis.range[1]"]
    check_if_update_necessary(start, end, display_start, display_end, old_contig, contig)
    if isinstance(GFF, pd.DataFrame):
        data = GFF
        data = data[(data["seqid"] == contig) & (data["start"] <= end) & (data["end"] >= start)]
        data["id"] = data.index
        data = data.to_dict('records')
    elif isinstance(GFF, str):
        db = gffutils.FeatureDB(GFF)
        data = db.region(start=start, end=end, seqid=contig)
        data = [{colname: getattr(row, colname) for colname in GFF_COLNAMES if colname != "attributes"} | {
            "attributes": str(row).split("\t")[-1], "id": row.id} for row in data]
    else:
        raise dash.exceptions.PreventUpdate

    return data


@callback(
    Output("coverage-graph", "figure", allow_duplicate=True),
    Input("trace-colors", "data"),
    State("coverage-graph", "figure"),

    prevent_initial_call=True
)
def my_callback(trace_colors, fig):
    if fig is None:
        raise dash.exceptions.PreventUpdate

    # Creating a Patch object
    patched_figure = dash.Patch()

    for idx, trace in enumerate(fig["data"]):
        if "legendgroup" in trace:
            color = trace_colors[trace["legendgroup"]]
            if trace["name"] == "Mean" or trace["name"] == "Median":
                patched_figure["data"][idx]["line"]["color"] = color
            else:
                a_color = css_color_to_rgba(color, 0.4)
                patched_figure["data"][idx]["fillcolor"] = a_color
    return patched_figure


@callback(
    Output("coverage-graph", "figure"),
    Output("internal-axis-range", "data"),
    Output("x-axis-range", 'data'),
    Output("current-contig", 'data'),
    Input("coverage-start", "value"),
    Input("coverage-end", "value"),
    Input("mode-switch", "value"),
    State("contig", "value"),
    State("current-contig", "data"),
    State("x-axis-range", "data"),
    State("autorange-y", "value"),
    State("coverage-graph", "figure"),
    State("trace-colors", "data"),

)
def update_coverage_plot(start, end, switch, contig, old_contig, axis_range, autorange_y, old_fig, trace_colors):
    design = COVERAGE_DESIGN
    coverage = COVERAGE_DATA
    display_start = axis_range["xaxis.range[0]"]
    display_end = axis_range["xaxis.range[1]"]
    print("updating figure", display_start, display_end, start, end, old_contig, contig)
    check_if_update_necessary(start, end, display_start, display_end, old_contig, contig)
    s = time.time()

    internal_window, step = _calc_internal_params(start, end)
    winsize = end - start
    print(internal_window, step)
    show_annotations = winsize <= 5000
    show_features = winsize <= 100000
    gff = GFF if isinstance(GFF, pd.DataFrame) else gffutils.FeatureDB(GFF)
    fig = plot_precomputed_coverage(
        design,
        contig=contig,
        coverages=coverage,
        wstart=internal_window[0],
        wend=internal_window[1],
        gff=gff,
        gff_name="gene_name",
        vertical_spacing=0,
        arrow_size=winsize * 0.01,
        step=step,
        show_annotations=show_annotations,
        show_features=show_features,
        type_colors=FEATURE_COLORS,
        colors=trace_colors,
        line_mapping=LINE_MAPPING
    )
    if not show_features:
        fig.add_annotation(
            text="Zoom in<br>to display<br>features", showarrow=False,
            row=2, col=1
        )
    if not switch:
        fig.update_layout(DARK_LAYOUT)
        linecolor = "white"
        fig.update_annotations(
            font=dict(color="black")
        )
    else:
        fig.update_layout(LAYOUT)
        linecolor = "black"

    fig.update_xaxes(
        showline=True,
        mirror=True,
        linecolor=linecolor,
        range=[start, end],
        minallowed=0,

    )
    fig.update_xaxes(
        row=3
    )
    fig.update_xaxes(
        zeroline=False,
        row=2,
    )
    fig.update_yaxes(
        showline=True,
        mirror=True,
        linecolor=linecolor

    )
    range_max = 4
    fig.update_yaxes(
        # showgrid=False,
        zeroline=False,
        range=[-1, range_max],
        row=2,
    )
    print(fig.layout["yaxis2"].maxallowed)
    if fig.layout["yaxis2"].maxallowed <= range_max:
        fig.update_yaxes(
            fixedrange=True,
            maxallowed=range_max,
            showgrid=False,
            showticklabels=False,
            row=2,
        )
        print(fig.layout["yaxis2"])

    fig.update_yaxes(
        minallowed=-10,
        row=1
    )
    fig.update_yaxes(
        minallowed=-10,
        row=3
    )
    if not autorange_y:
        for axis in ["", "3"]:
            axis = f"yaxis{axis}"
            fig["layout"][axis]["range"] = old_fig["layout"][axis]["range"]

    fig.update_layout(dragmode="pan", uirevision=True)
    fig.update_shapes(line=dict(color=linecolor))
    lout = {'xaxis.range[0]': start, 'xaxis.range[1]': end, }
    e = time.time()
    print(f"no honestly im updating {e - s}")

    return fig, internal_window, lout, contig


layout = get_layout()
