
import os
import dash
from dash import Dash, html, dcc, clientside_callback, Input, Output, State, ALL, ClientsideFunction
import dash_bootstrap_components as dbc
import DEplots.dashboard

FILEDIR = os.path.dirname(os.path.abspath(__file__))
ASSETS_DIR = os.path.join(FILEDIR, "assets")
assert os.path.exists(ASSETS_DIR)

CONFIG = DEplots.dashboard.CONFIG

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
                dcc.Location(id='url', refresh="callback-nav"),

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
                        dbc.NavItem(dbc.NavLink(f"{page['name']}", href=page["relative_path"], id={"type": f"nav-item", "index": idx}), className="p-1") for
                        idx, page in enumerate(dash.page_registry.values())
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
        className="ufr-navbar shadow w-100", style={"position": "fixed", "z-index": "99999"}
    )
    return navbar


def _get_footer():
    div = html.Footer(
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
                                        html.Li(html.A(html.I(className="fa-solid fa-2xl fa-envelope"), target="_blank", href=f"mailto:{CONFIG['email']}", className="text-light")),
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
    return div

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
            ),
            _get_footer()
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


@app.callback(
    Output({'type': 'nav-item', 'index': ALL}, 'active'),
    Input("url", "pathname"),
    State({'type': 'nav-item', 'index': ALL}, 'href'),
)
def update_active_nav(url, state):
    d = [url == href for href in state]
    return d


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

