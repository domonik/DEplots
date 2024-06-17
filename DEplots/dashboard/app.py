
import os
import dash
from dash import Dash, html, dcc, clientside_callback, Input, Output, State
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


color_mode_switch = html.Span(
    [
        dbc.Label(className="fa fa-xl fa-moon", html_for="switch", style={"vertical-align": "0 !important"}),
        dbc.Switch( id="mode-switch", value=True, className="d-inline-block ms-3 fs-4", persistence=True),
        dbc.Label(className="fa fa-xl fa-sun", html_for="switch", style={"vertical-align": "0 !important"}),
    ]
)




def get_navbar():
    navbar = dbc.Navbar(
        dbc.Container(
            [
                html.A(
                    dbc.Row(
                        [

                            dbc.Col(html.Img(
                                src=app.get_asset_url("BioInfLogo.png"),
                                height="50px"),
                            ),
                            dbc.Col(dbc.NavbarBrand("DESeq Explorer", className="ms-2"), className="d-flex align-items-center"),
                        ]

                    )

                ),
                dbc.NavbarToggler(id="navbar-toggler", n_clicks=0),

                dbc.Collapse(
                    [
                        dbc.NavItem(dbc.NavLink(f"{page['name']}", href=page["relative_path"]), className="p-1") for
                        page in dash.page_registry.values()
                    ],
                    is_open=False,
                    navbar=True,
                    id="navbar-collapse2",

                ),
                dbc.Collapse(
                    [color_mode_switch] + [
                        dbc.NavItem(
                            html.Img(
                                src="https://cd.uni-freiburg.de/wp-content/uploads/2022/09/ufr-logo-white-2.png",
                                height="50px"),
                            className="ps-0 ps-md-3"

                        )
                    ],
                    id="navbar-collapse",
                    className="justify-content-end",
                    is_open=False,
                    navbar=True,
                ),

                # dbc.Row(
                #     [
                #
                #         dbc.Col(
                #             dbc.NavbarToggler(id="navbar-toggler", n_clicks=0),
                #             className="d-flex d-md-none"
                #         ),
                #         dbc.Col(
                #             [
                #                 dbc.NavItem(dbc.NavLink(f"{page['name']}", href=page["relative_path"]), className="p-1") for page in
                #                 dash.page_registry.values()
                #             ],
                #             width=3, md=3, className="d-flex align-items-center"
                #         ),
                #         dbc.Col(html.H2("DESeq Explorer", className="text-md-center"), width=3, align="center",),
                #
                #         dbc.Col(
                #             color_mode_switch,
                #             width=4, className="d-flex align-items-center"
                #         )
                #     ],
                #     #className="w-100 ",
                #
                # ),



            ],
            fluid=True,
            className="dbc text-light"

        ),
        dark=True,
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


@app.callback(
    Output("navbar-collapse", "is_open"),
    Output("navbar-collapse2", "is_open"),
    [Input("navbar-toggler", "n_clicks")],
    [
        State("navbar-collapse", "is_open"),
    ],
)
def toggle_navbar_collapse(n, is_open):
    if n:
        return not is_open, not is_open
    return is_open, is_open

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

