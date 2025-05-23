import dash
from dash import html, callback, Input, Output
from DEplots.dashboard import CONFIG
from dash import dcc
import os
import re
import base64
import dash_bootstrap_components as dbc

dash.register_page(__name__, path='/')

if CONFIG["welcome_files"]:
    WELCOME_FILES = CONFIG["welcome_files"]

else:
    file = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "defaultWelcome.md")
    WELCOME_FILES = [file]


def replace_relative_links(markdown_text, markdown_path, python_path):
    # Extract directory path of the markdown file
    markdown_dir = os.path.dirname(markdown_path)

    # Extract directory path of the python file
    python_dir = os.path.dirname(python_path)

    # Replace relative links in img tags
    def replace_link(match):
        link = match.group(1)
        # Make the link path relative to the Python file
        img_path = os.path.join(markdown_dir, link)

        # Determine the file extension
        ext = os.path.splitext(img_path)[1].lower()

        # Read image file
        with open(img_path, 'rb') as f:
            img_data = f.read()

        # Encode based on image format
        if ext == '.svg':
            img_data_encoded = base64.b64encode(img_data).decode()
            img_src = f'data:image/svg+xml;base64,{img_data_encoded}'
        elif ext in ('.png', '.jpg', '.jpeg'):
            img_data_encoded = base64.b64encode(img_data).decode()
            img_src = f'data:image/{ext[1:]};base64,{img_data_encoded}'
        else:
            raise ValueError("Unsupported image format")

        return f'src="{img_src}"'

    # Regular expression to match HTML img tags and extract src attribute
    img_regex = r'src="([^"]+)"'

    # Replace links in the markdown text
    new_markdown_text = re.sub(img_regex, replace_link, markdown_text)

    return new_markdown_text


WELCOME_TEXTS = []

for file in WELCOME_FILES:
    with open(file) as handle:
        text = handle.read()
    md_file = os.path.abspath(file)
    assert os.path.exists(md_file)
    text = replace_relative_links(text, md_file, os.path.abspath(__file__))
    WELCOME_TEXTS.append(text)


def welcome_layout():
    welcome = [

        dbc.Col(
            [
                dbc.Card(
                    [
                        html.Div(
                            [
                                dcc.Markdown(text, dangerously_allow_html=True, ),
                            ]
                        )
                    ]
                    , className="shadow p-2", style={"font-size": "20px"})

            ],
            width=12, lg=6, className="py-1"
        ) for text in WELCOME_TEXTS
    ]

    return welcome



layout = html.Div(
        welcome_layout(),
    style={
        'minHeight': '85vh',  # Ensures the div takes at least the full viewport height
        'width': '100%',  # Full width
        'padding': '20px'  # Optional: adds some spacing
    },
    className = "row p-1 justify-content-between" if len(WELCOME_TEXTS) > 1 else "row p-1 justify-content-around"


)
