import pandas as pd
import os
import json
import mmap
import numpy as np
import plotly
from plotly import graph_objs as go
from scipy.stats import linregress
from find_top_genes import sort_genes
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
from textwrap import dedent as d


def generate_fig(value_dict):
    trace = go.Scatter(
        x=value_dict["conc_nonzero"],
        y=value_dict["reads_log"],
        mode='markers',
        showlegend=False,
        marker=dict(
            color='#1f77b4'
        ),
        xaxis='x1',
        yaxis='y1'
    )
    trace2 = go.Scatter(
        x=value_dict["fit_vals_x"],
        y=value_dict["fit_vals_y"],
        mode='lines',
        showlegend=False,
        xaxis='x1',
        yaxis='y1'
    )
    trace3 = go.Scatter(
        x=value_dict["conc_nonzero"],
        y=value_dict["reads_norm_log"],
        mode='markers',
        showlegend=False,
        marker=dict(
            color='#1f77b4'
        ),
        xaxis='x2',
        yaxis='y2'
    )
    trace4 = go.Scatter(
        x=value_dict["fit_vals_norm_x"],
        y=value_dict["fit_vals_norm_y"],
        mode='lines',
        showlegend=False,
        xaxis='x2',
        yaxis='y2'
    )
    layout = go.Layout(
        annotations=[
            dict(
                x=-3,
                y=np.max(value_dict["reads_log"]),
                xref='x1',
                yref='y1',
                text=f'R-squared: {value_dict["r_value"]**2:0.2f}<br>Stderr: {value_dict["std_err"]:02f}',
                showarrow=False,
                font=dict(
                    size=14
                )
            ),
            dict(
                x=-3,
                y=np.max(value_dict["reads_norm_log"]),
                xref='x2',
                yref='y2',
                text=f'R-squared: {value_dict["r_value_norm"]**2:0.2f}<br>Stderr: {value_dict["std_err_norm"]:02f}',
                showarrow=False,
                font=dict(
                    size=14
                )
            )
        ],
        xaxis=dict(
            title='log(Concentration)',
            domain=[0, 0.45],
            anchor='y1'
        ),
        xaxis2=dict(
            title='log(Concentration)',
            domain=[0.55, 1],
            anchor='y2'
        ),
        yaxis=dict(
            title='log(Reads)',
            anchor='x1'
        ),
        yaxis2=dict(
            title='log(IC Normalized Reads)',
            anchor='x2'
        )
    )
    fig = go.Figure([trace, trace2, trace3, trace4], layout)
    return fig


def generate_ctrls_fig(value_dict):
    trace = go.Scatter(
        x=value_dict["conc_nonzero"],
        y=value_dict["ctrls_nonzero"][:, 0],
        mode='markers',
        showlegend=False,
        marker=dict(
            color='#1f77b4'
        )
    )
    trace1 = go.Scatter(
        x=value_dict["conc_nonzero"],
        y=value_dict["ctrls_nonzero"][:, 1],
        mode='markers',
        showlegend=False,
        marker=dict(
            color='#1f77b4'
        )
    )
    fig = go.Figure([trace, trace1])
    return fig


def filter_data(count_json, gene, organism, fit_range, fit_range_norm):
    org_json = count_json[str(organism)]
    out_dict = {}
    conc = np.array(org_json["Gene Counts"][gene]["Concentration"])
    # D1 = 7.36E6 cfu/ml for Kpneumoniae
    # D1 = 1.16E7 cfu/ml for Bpertussis
    # D1 = 1.88E6 cfu/ml for Saureus
    reads, ctrls = org_json["Gene Counts"][gene]['Read Counts'], org_json['Ctrl Counts']
    reads_nonzero_idx = np.argwhere(np.array(reads) > 0).reshape(-1)
    reads_nonzero = np.array(reads)[reads_nonzero_idx]
    reads_nonzero[reads_nonzero == 0] = 1
    ctrls_nonzero = np.array(ctrls)[reads_nonzero_idx]
    conc_nonzero = conc[reads_nonzero_idx]
    reads_log = np.log10(reads_nonzero)
    fit_idx = np.argwhere((conc_nonzero >= fit_range[0]) & (conc_nonzero <= fit_range[1])).reshape(-1)
    slope, intercept, r_value, p_value, std_err = linregress(conc_nonzero[fit_idx], reads_log[fit_idx])
    fit_vals_y = slope * conc_nonzero[fit_idx] + intercept
    fit_vals_x = conc_nonzero[fit_idx]
    # Normalized
    reads_norm = reads_nonzero / ctrls_nonzero
    reads_norm_log = np.log10(reads_norm)
    # reads_norm_log_mean = np.mean(reads_norm_log, axis=1)
    fit_norm_idx = np.argwhere((conc_nonzero >= fit_range_norm[0]) & (conc_nonzero <= fit_range_norm[1])).reshape(-1)
    slope_norm, intercept_norm, r_value_norm, p_value_norm, std_err_norm = \
        linregress(conc_nonzero[fit_norm_idx], reads_norm_log[fit_norm_idx])
    fit_vals_norm_y = slope_norm * conc_nonzero[fit_norm_idx] + intercept_norm
    fit_vals_norm_x = conc_nonzero[fit_norm_idx]
    out_dict.update({
        "conc_nonzero": conc_nonzero,
        "reads_log": reads_log,
        "reads_norm_log": reads_norm_log,
        "fit_vals_y": fit_vals_y,
        "fit_vals_x": fit_vals_x,
        "fit_vals_norm_y": fit_vals_norm_y,
        "fit_vals_norm_x": fit_vals_norm_x,
        "r_value": r_value,
        "r_value_norm": r_value_norm,
        "std_err": std_err,
        "std_err_norm": std_err_norm,
        "ctrls_nonzero": ctrls_nonzero
    })
    return out_dict


orgs_taxid = {
    # 'Hinfluenzae': 727,
    'Bpertussis': 520,
    'Kpneumoniae': 573,
    'Saureus': 1280,
    'Nfarcinica': 37329,
    # 'Nwallacei': 480035
}

with open('./data/lod_counts_dna_raw.json') as dna_file:
    dna_counts = json.load(dna_file)

with open('./data/lod_counts_rna_raw.json') as rna_file:
    rna_counts = json.load(rna_file)

app = dash.Dash(__name__)

app.layout = html.Div([
    html.Div([
        dcc.Markdown(d("""
            ### LOD Linearity
        """))
    ], className='row'),
    html.Div([
        html.Span([
            html.P([
                'Organism'
            ]),
            dcc.Dropdown(
                id='organism-select',
                options=[{'label': key, 'value': val} for (key, val) in orgs_taxid.items()],
                value=573
            )
        ], className='two columns'),
        html.Div([
            html.P([
                'Library Type'
            ]),
            dcc.Dropdown(
                id='libtype-select',
                options=[{'label': val, 'value': val} for val in ['DNA', 'RNA']],
                value='RNA'
            )
        ], className='two columns'),
        html.Div([
            html.P([
                'Gene'
            ]),
            dcc.Dropdown(
                id='gene-select',
                options=[{'label': '0', 'value': '0'}],
                value='0'
            )
        ], className='two columns')
    ], className='row'),
    html.Div([
        dcc.Graph(
            id='graph'
        )
    ], className='row'),
    html.Div([
        html.Div([
            dcc.RangeSlider(
                id='fit_range',
                marks={-5: '-5', -4: '-4', -3: '-3', -2: '-2', -1: '-1', 0: '0'},
                min=-5,
                max=0,
                value=[-5, 0]
            )
        ], className='six columns', style={'padding': '25px 80px 50px 120px'}),
        html.Div([
            dcc.RangeSlider(
                id='fit_range_norm',
                marks={-5: '-5', -4: '-4', -3: '-3', -2: '-2', -1: '-1', 0: '0'},
                min=-5,
                max=0,
                value=[-5, 0]
            )
        ], className='six columns', style={'padding': '25px 120px 50px 80px'})
    ], className='row')
], style={'padding': '10px 50px 50px 50px'})


@app.callback(
    Output('graph', 'figure'),
    [Input('organism-select', 'value'),
     Input('libtype-select', 'value'),
     Input('gene-select', 'value'),
     Input('fit_range', 'value'),
     Input('fit_range_norm', 'value')]
)
def update_graph(orgSelect, libSelect, geneSelect, fit_range, fit_range_norm):
    if libSelect == 'DNA':
        count_data = dna_counts
    elif libSelect == 'RNA':
        count_data = rna_counts
    value_dict = filter_data(count_data, geneSelect, orgSelect, fit_range, fit_range_norm)
    figure = generate_fig(value_dict)
    return figure


@app.callback(
    [Output('gene-select', 'value'),
     Output('gene-select', 'options')],
    [Input('organism-select', 'value'),
     Input('libtype-select', 'value')]
)
def update_gene_select(orgSelect, libSelect):
    if libSelect == 'DNA':
        count_data = dna_counts
        keys_list = list(count_data[str(orgSelect)]["Gene Counts"].keys())
        genes_sorted = sort_genes(count_data, orgSelect)
        options =  [{'label': val, 'value': val} for val in genes_sorted]
        value = keys_list[0]
    elif libSelect == 'RNA':
        count_data = rna_counts
        options = [{'label': '0', 'value': '0'}]
        value='0'
    return value, options


if __name__ == '__main__':
    app.run_server(debug=True)
