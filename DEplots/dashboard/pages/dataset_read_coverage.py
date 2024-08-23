import dash
import pandas as pd
from dash import Input, Output, State, html, dcc, callback, Patch
import dash_bootstrap_components as dbc
from dash.dash_table.Format import Format, Scheme
from pandas.api.types import is_numeric_dtype
from DEplots.runComparison import plot_gene_among_conditions, upset_plot_from_deseq
from DEplots.dashboard import DEFAULT_PLOTLY_COLORS, DEFAULT_PLOTLY_COLORS_LIST, LAYOUT, DARK_LAYOUT, UP_COLOR_LIGHT, \
    UP_COLOR_DARK, DOWN_COLOR_LIGHT, DOWN_COLOR_DARK, COVERAGE_DATA, COVERAGE_DESIGN, GFF
import numpy as np
from DEplots.readCount import plot_precomputed_coverage


if COVERAGE_DATA is not None:
    dash.register_page(__name__, path='/coverage', name="Read Coverage")

SETTINGS_ROW = "px-2"
SETTINGS_COL = "py-2 d-flex align-items-center"


def coverage_settings_card():
    contigs = list(COVERAGE_DATA.keys())
    coverage_plot_card = dbc.Col(
        dbc.Card(
            [
                dbc.CardHeader(
                    dbc.Row(
                        [
                            dbc.Col(html.H5("Settings"), width=6),

                        ]

                    ),

                ),
                dbc.Col(
                    [
                        dbc.Row(
                            [
                                dbc.Col(html.Span("Sequence"), width=6, md=3, className=SETTINGS_COL),
                                dbc.Col(dcc.Dropdown(id="contig", options=contigs, value=contigs[0], clearable=False, style={"width": "100%"}),
                                        width=6, className=SETTINGS_COL,
                                        md=3),
                                dbc.Col(html.Span("Autorange Y"), width=6, md=3, className=SETTINGS_COL),
                                dbc.Col(dbc.Switch(id="autorange-y", value=True, className="fs-5"), width=6, md=3,
                                        className="d-flex justify-content-end " + SETTINGS_COL),

                            ],
                            className=SETTINGS_ROW

                        ),

                        dbc.Row(
                            [
                                dbc.Col(html.Span("Start"), width=6, md=3, className=SETTINGS_COL),
                                dbc.Col(dbc.Input(id="coverage-start", type="number", value=0, min=0, step=1, className=SETTINGS_COL),
                                        width=6, md=3),
                                dbc.Col(html.Span("End"), width=6, md=3,className=SETTINGS_COL),
                                dbc.Col(dbc.Input(id="coverage-end", type="number", value=2000, min=0, step=1, className=SETTINGS_COL),
                                        width=6, md=3),
                            ],
                            className=SETTINGS_ROW

                        ),

                    ],
                    className="p-2"

                ),



            ],
            className="shadow",
        ),
        width=12, md=6
    )

    return coverage_plot_card


def coverage_card():
    coverage_plot_card = dbc.Col(
        dbc.Card(
            [
                dbc.CardHeader(
                    dbc.Row(
                        [
                            dbc.Col(html.H5("Coverage"), width=6),
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
                dbc.Col(
                    dcc.Graph(
                        id="coverage-graph",
                        style={"resize": "vertical", "overflow": "auto"},
                        config={
                            "scrollZoom": True,

                        }
                    ), width=12,

                ),

            ],
            className="shadow",
        ),
        width=12
    )

    return coverage_plot_card


def get_layout():
    lout = html.Div(
        [
            dcc.Store(id="x-axis-range", data={"xaxis.range[0]": 1, "xaxis.range[1]": 2000, "yaxis2.range[0]": -0.5, "y2axis.range[1]": 5.5}),
            dcc.Store(id="internal-axis-range", data=[0, 2000]),
            dbc.Container(
                [
                    dbc.Row(
                        coverage_card(),
                        className="py-1"

                    ),
                    dbc.Row(
                        coverage_settings_card(),
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
    if np.isclose(old_winsize, new_winsize): # only shifts window
        return False
    else:
        return True




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

        is_zoom = _is_zoom_relayout(old_lout, new_lout)
        if is_zoom:
            print("zoomed")
            print(old_lout)
            print(new_lout)
            new_lout["xaxis.range[0]"] = -1  # This forces the figure to replot

            return new_lout, start, end
        else: # now we need to determine if we moved outside the window_range
            if new_lout["xaxis.range[1]"] > internal_xaxis[1] or new_lout["xaxis.range[0]"] < internal_xaxis[0]:
                print("scrolled out of window")
                new_lout["xaxis.range[0]"] = -1 # This forces the figure to replot
                print(new_lout)
                return new_lout,  start, end
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
    Output("coverage-start", "value"),
    Output("coverage-end", "value"),
    Input("contig", "value")
)
def update_via_contig(contig):
    start=0
    end = min(COVERAGE_DATA[contig]["+"].shape[-1], 2000)
    return start, end


@callback(
    Output("coverage-graph", "figure"),
    Output("internal-axis-range", "data"),
    Output("x-axis-range", 'data'),
    Input("coverage-start", "value"),
    Input("coverage-end", "value"),
    Input("mode-switch", "value"),
    State("contig", "value"),
    State("x-axis-range", "data"),
    State("autorange-y", "value"),
    State("coverage-graph", "figure"),
)
def update_coverage_plot(start, end, switch, contig, axis_range, autorange_y, old_fig):
    design = COVERAGE_DESIGN
    coverage = COVERAGE_DATA
    display_start = axis_range["xaxis.range[0]"]
    display_end = axis_range["xaxis.range[1]"]
    print(display_start, display_end, start, end)
    if dash.ctx.triggered_id != "mode-switch":
        if int(display_start) == int(start) and int(display_end) == int(end):
            raise dash.exceptions.PreventUpdate
        else:
            display_start = start
            display_end = end
        if display_start is None or display_end is None or display_start >= display_end:
            raise dash.exceptions.PreventUpdate



    internal_window, step = _calc_internal_params(display_start, display_end)
    winsize = display_end - display_start
    print(internal_window, step)
    show_annotations = winsize <= 5000
    fig = plot_precomputed_coverage(
        design,
        contig=contig,
        coverages=coverage,
        wstart=internal_window[0],
        wend=internal_window[1],
        gff=GFF,
        gff_name="locus_tag",
        vertical_spacing=0,
        arrow_size=winsize * 0.01,
        step=step,
        show_annotations=show_annotations
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
        range=[display_start, display_end],
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
    range_max = 5.5
    fig.update_yaxes(
        showgrid=False,
        zeroline=False,
        range=[-0.5, range_max],
        row=2,
    )
    if fig.layout["yaxis2"].maxallowed <= range_max:
        fig.update_yaxes(
            fixedrange = True,
            row=2
        )
    fig.update_yaxes(
        minallowed=-10,
        row=1
    )
    fig.update_yaxes(
        minallowed=-10,
        row=2
    )
    if not autorange_y:
        for axis in ["", "2"]:
            axis = f"yaxis{axis}"
            fig["layout"][axis]["range"] = old_fig["layout"][axis]["range"]

    fig.update_layout(dragmode="pan", uirevision=True)
    fig.update_shapes(line=dict(color=linecolor))
    lout = {'xaxis.range[0]': display_start, 'xaxis.range[1]': display_end,}
    return fig, internal_window, lout




layout = get_layout()