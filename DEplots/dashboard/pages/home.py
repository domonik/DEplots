import dash
from dash import html, callback, Input, Output

dash.register_page(__name__, path='/')

layout = html.Div([
    html.H1('This is our Home page', ),
    html.Div('This is our Home page content.', id="header"),
])

@callback(
    Output("header", "children"),
    Input("mode-switch", "value")
)
def change(mode):
    return "Dark Homepage" if mode else "Light Homepage"
