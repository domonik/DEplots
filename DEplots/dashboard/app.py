
import os
import dash
from dash import Dash, html, dcc, clientside_callback, Input, Output
import dash_bootstrap_components as dbc

FILEDIR = os.path.dirname(os.path.abspath(__file__))
ASSETS_DIR = os.path.join(FILEDIR, "assets")
assert os.path.exists(ASSETS_DIR)

dbc_css = "https://cdn.jsdelivr.net/gh/AnnMarieW/dash-bootstrap-templates/dbc.min.css"
print("Defining APP")

app = Dash(
    __name__,
    external_stylesheets=["custom.css", dbc.icons.FONT_AWESOME, dbc_css],
    assets_folder=ASSETS_DIR,
    prevent_initial_callbacks='initial_duplicate',
    use_pages=True
)
print("APP defined")


color_mode_switch = dbc.Row(
    [
        dbc.Col(html.I(className="fa fa-xl fa-moon", ), className="d-flex justify-content-end align-items-center", width=3, md=1),
        dbc.Col(dbc.Switch(id="mode-switch", className="ms-3 fs-4", persistence=True), className="d-flex justify-content-center align-items-center", width=3, md=1),
        dbc.Col(html.I(className="fa fa-xl fa-sun"), className="d-flex justify-content-start align-items-center", width=3, md=1 ),
    ], justify="end", className="w-100"
)




def get_navbar():
    navbar = dbc.Navbar(
        dbc.Container(
            [
                dbc.Row(
                    [
                        dbc.Col(
                            dbc.Row(
                                [
                                    dbc.Col(html.Img(src="https://cd.uni-freiburg.de/wp-content/uploads/2022/09/ufr-logo-white-2.png", height="50px")),

                                ] + [
                                    dbc.Col(dcc.Link(f"{page['name']}", href=page["relative_path"])) for page in dash.page_registry.values()
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
    return navbar


def get_layout():
    layout = html.Div(
        [
            get_navbar(),
            dcc.Store(data={}, id="volcano-highlight-ids"),
            dcc.Store(id="enrich-term"),
            dbc.Container(
                dash.page_container,
                fluid=True,
                style={"padding-top": "65px"}
            )
        ],
    )
    return layout


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

