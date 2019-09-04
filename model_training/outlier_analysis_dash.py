import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_daq as daq
from dash.dependencies import Input, Output
from textwrap import dedent


app = dash.Dash(__name__)
app.layout = html.Div([
    html.Div([
        dcc.Markdown(dedent("""
            ### Organism Coverage
        """))
    ], className='row'),
    html.Div([
        html.Div([
            html.P([
                'Organism'
            ]),
            dcc.Dropdown(
                id='org-select',
                options=[{'label': label, 'value': value} for label, value in
                    zip(list(org_cov.keys()), [ncbi.get_name(key) for key in org_cov.keys()])],
                value=287
            )
        ], className='two columns'),
        html.Div([
            html.P([
                '10<sup>1</sup>'
            ]),
            dcc.Graph(id='graph10-1')
        ], className='row'),
        html.Div([
            html.P([
                '10<sup>1.5</sup>'
            ]),
            dcc.Graph(id='graph10-15')
        ], className='row'),
        html.Div([
            html.P([
                '10<sup>2</sup>'
            ]),
            dcc.Graph(id='graph10-2')
        ], className='row'),
        html.Div([
            html.P([
                '10<sup>2.5</sup>'
            ]),
            dcc.Graph(id='graph10-25')
        ], className='row'),
        html.Div([
            html.P([
                '10<sup>3</sup>'
            ]),
            dcc.Graph(id='graph10-3')
        ], className='row'),
        html.Div([
            html.P([
                '10<sup>3.5</sup>'
            ]),
            dcc.Graph(id='graph10-35')
        ], className='row'),
        html.Div([
            html.P([
                '10<sup>4</sup>'
            ]),
            dcc.Graph(id='graph10-4')
        ], className='row'),
    ])
])