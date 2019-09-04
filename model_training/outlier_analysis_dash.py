import json
import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_daq as daq
from dash.dependencies import Input, Output
from textwrap import dedent
from ncbi_taxonomy_utils import ncbi_taxonomy
import numpy as np
import plotly.graph_objects as go

ncbi = ncbi_taxonomy()
with open('model_training/coverage_strings.json') as file:
    org_cov = json.load(file)
colors = [
    '#1f77b4',
    '#ff7f0e',
    '#2ca02c',
    '#d62728',
    '#9467bd',
    '#8c564b',
    '#e377c2',
    '#7f7f7f',
    '#bcbd22',
    '#17becf'
]

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
                    zip([ncbi.get_name(int(key)) for key in org_cov.keys()], list(org_cov.keys()))],
                value="287"
            )
        ], className='row'),
        html.Div([
            # html.P([
            #     '10^1'
            # ]),
            dcc.Graph(id='graph10-1',
                config={
                    'displayModeBar': False
            })
        ], className='row'),
        html.Div([
            # html.P([
            #     '10^1.5'
            # ]),
            dcc.Graph(id='graph10-15',
                config={
                    'displayModeBar': False
            })
        ], className='row'),
        html.Div([
            # html.P([
            #     '10^2'
            # ]),
            dcc.Graph(id='graph10-2',
                config={
                    'displayModeBar': False
            })
        ], className='row'),
        html.Div([
            # html.P([
            #     '10^2.5'
            # ]),
            dcc.Graph(id='graph10-25',
                config={
                    'displayModeBar': False
            })
        ], className='row'),
        html.Div([
            # html.P([
            #     '10^3'
            # ]),
            dcc.Graph(id='graph10-3',
                config={
                    'displayModeBar': False
            })
        ], className='row'),
        html.Div([
            # html.P([
            #     '10^3.5'
            # ]),
            dcc.Graph(id='graph10-35',
                config={
                    'displayModeBar': False
            })
        ], className='row'),
        html.Div([
            # html.P([
            #     '10^4'
            # ]),
            dcc.Graph(id='graph10-4',
                config={
                    'displayModeBar': False
            })
        ], className='row'),
    ])
])


@app.callback(
    [
        Output('graph10-1', 'figure'),
        Output('graph10-15', 'figure'),
        Output('graph10-2', 'figure'),
        Output('graph10-25', 'figure'),
        Output('graph10-3', 'figure'),
        Output('graph10-35', 'figure'),
        Output('graph10-4', 'figure'),
    ],
    [
        Input('org-select', 'value')
    ]
)
def update_graphs(orgSelect):
    cov_dict = org_cov[orgSelect]
    dilutions = list(set(cov_dict['Dilution']))
    dilutions.sort()
    graph_ls = ['graph10-1', 'graph10-15', 'graph10-2', 'graph10-25',
                'graph10-3', 'graph10-35', 'graph10-4']
    graph_dict = {}
    for idx, (dil, graph) in enumerate(zip(dilutions, graph_ls)):
        dil_idx = np.argwhere(np.array(cov_dict['Dilution']) == dil)
        coverage_ls = []
        for i in dil_idx.flatten():
            coverage_ls.append([float(el) for el in cov_dict['Coverage'][i]])
        data = []
        for coverage in coverage_ls:
            data.append(
                go.Scatter(
                    y=coverage,
                    mode='lines',
                    line=dict(
                        color=colors[idx]
                    ),
                    showlegend=False
                )
            )
        layout = go.Layout(
            height=100,
            margin=go.layout.Margin(
                l=100,
                r=0,
                b=0,
                t=0,
                pad=0
            )
        )
        fig=go.Figure(data, layout)
        graph_dict.update({graph: fig})
    return graph_dict[graph_ls[0]], graph_dict[graph_ls[1]], graph_dict[graph_ls[2]], graph_dict[graph_ls[3]], graph_dict[graph_ls[4]], graph_dict[graph_ls[5]], graph_dict[graph_ls[6]]


if __name__ == '__main__':
    app.run_server(debug=True, port=9080)