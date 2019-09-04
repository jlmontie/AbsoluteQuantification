import pandas as pd
import os
import json
import mmap
import numpy as np
import plotly
from plotly import graph_objs as go
from scipy.stats import linregress, mode
from find_top_genes import sort_genes
import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_daq as daq
from dash.dependencies import Input, Output
from textwrap import dedent


def generate_fig(value_dict):
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
    data = []
    annotations = []
    slope_ls = []
    slope_norm_ls = []
    r2_ls = []
    r2_norm_ls = []
    for idx, org in enumerate(value_dict.keys()):
        org_dict = value_dict[org]
        slope_ls.append(org_dict["slope"])
        slope_norm_ls.append(org_dict["slope_norm"])
        r2_ls.append(org_dict["r_value"])
        r2_norm_ls.append(org_dict["r_value_norm"])
        trace = go.Scatter(
            x=org_dict["conc_nonzero"],
            y=org_dict["reads_log"],
            mode='markers',
            name=taxid_orgs[int(org)],
            marker=dict(
                color=colors[idx]
            ),
            xaxis='x1',
            yaxis='y1',
            legendgroup=idx
        )
        trace2 = go.Scatter(
            x=org_dict["fit_vals_x"],
            y=org_dict["fit_vals_y"],
            mode='lines',
            line=dict(
                color=colors[idx]
            ),
            showlegend=False,
            name=taxid_orgs[int(org)],
            xaxis='x1',
            yaxis='y1',
            legendgroup=idx
        )
        trace3 = go.Scatter(
            x=org_dict["conc_nonzero"],
            y=org_dict["reads_norm_log"],
            mode='markers',
            showlegend=False,
            name=taxid_orgs[int(org)],
            marker=dict(
                color=colors[idx]
            ),
            xaxis='x2',
            yaxis='y2',
            legendgroup=idx
        )
        trace4 = go.Scatter(
            x=org_dict["fit_vals_norm_x"],
            y=org_dict["fit_vals_norm_y"],
            mode='lines',
            line=dict(
                color=colors[idx]
            ),
            showlegend=False,
            name=taxid_orgs[int(org)],
            xaxis='x2',
            yaxis='y2',
            legendgroup=idx
        )
        data.extend([trace, trace2, trace3, trace4])
    annotations.extend([
        dict(
            x=-2,
            y=np.max([np.max(value_dict[key]['reads_log']) for key in value_dict]),
            xref='x1',
            yref='y1',
            text=f'Mean slope: {np.mean(slope_ls):0.2f}+/-{np.max(slope_ls) - np.min(slope_ls) / 2:0.2f}<br>Mean R<sup>2</sup>: {np.mean(r2_ls)**2:0.2f}',
            showarrow=False,
            font=dict(
                size=14,
            )
        ),
        dict(
            x=4,
            y=np.max([np.max(value_dict[key]['reads_norm_log']) for key in value_dict]),
            xref='x2',
            yref='y2',
            text=f'Mean slope: {np.mean(slope_norm_ls):0.2f} +/- {(np.max(slope_norm_ls) - np.min(slope_norm_ls)) / 2:0.2f}<br>Mean R<sup>2</sup>: {np.median(r2_norm_ls)**2:0.2f}',
            showarrow=False,
            font=dict(
                size=14,
            )
        )
    ])
    layout = go.Layout(
        annotations = annotations,
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
            title='log(Coverage)',
            anchor='x1'
        ),
        yaxis2=dict(
            title='log(IC Normalized Coverage)',
            anchor='x2'
        )
    )
    fig = go.Figure(data, layout)
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


def filter_data(count_dict, copyNumber):#, fit_range, fit_range_norm):
    out_dict = {}
    for org in count_dict.keys():
        org_dict = count_dict[str(org)]
        conc = np.array(org_dict["Concentration"])
        # D1 = 2.9E7 cfu/ml for Bpertussis
        # D1 = 1.84E7 cfu/ml for Kpneumoniae
        # D1 = 4.7E6 cfu/ml for Saureus
        # D1 = 8.0E4 cfu/ml for Nfarcinica
        reads, ctrls = org_dict['Read Counts'], org_dict['Ctrl Counts']
        print(reads)
        if org in copy_numbers and copyNumber:
            copy_number_normalizer = copy_numbers[org]['copies']
        else:
            copy_number_normalizer = 1
        reads = np.array(reads) / copy_number_normalizer
        reads = np.array(reads)
        reads_nonzero_idx = np.argwhere((reads > 0)).reshape(-1)
        reads_nonzero = reads[reads_nonzero_idx]

        print(reads_nonzero)
        reads_nonzero[reads_nonzero == 0] = 1
        ctrls_nonzero = np.array(ctrls)[reads_nonzero_idx]
        conc_nonzero = conc[reads_nonzero_idx]
        reads_log = np.log10(reads_nonzero)
        slope, intercept, r_value, p_value, std_err = linregress(conc_nonzero, reads_log)
        fit_vals_y = slope * conc_nonzero + intercept
        fit_vals_x = conc_nonzero
        # Normalized
        reads_norm = reads_nonzero / ctrls_nonzero
        reads_norm_log = np.log10(reads_norm)

        slope_norm, intercept_norm, r_value_norm, p_value_norm, std_err_norm = \
            linregress(conc_nonzero, reads_norm_log)
        fit_vals_norm_y = slope_norm * conc_nonzero + intercept_norm
        fit_vals_norm_x = conc_nonzero
        org_dict.update({
            "conc_nonzero": conc_nonzero,
            "reads_log": reads_log,
            "reads_norm_log": reads_norm_log,
            "slope": slope,
            "fit_vals_y": fit_vals_y,
            "fit_vals_x": fit_vals_x,
            "slope_norm": slope_norm,
            "fit_vals_norm_y": fit_vals_norm_y,
            "fit_vals_norm_x": fit_vals_norm_x,
            "r_value": r_value,
            "r_value_norm": r_value_norm,
            "std_err": std_err,
            "std_err_norm": std_err_norm,
            "ctrls_nonzero": ctrls_nonzero
        })
        out_dict.update({org: org_dict})
    return out_dict


with open('./data/community_standard_panel.json') as panel_file:
    taxid_orgs = json.load(panel_file)
taxid_orgs = {int(k):v for k,v in taxid_orgs.items()}

with open('data/community_std_counts/community_std_counts_16s_a.json') as a_file:
    a = json.load(a_file)

with open('data/community_std_counts/community_std_counts_16s_b.json') as b_file:
    b = json.load(b_file)

with open('data/community_std_counts/community_std_counts_16s_c.json') as c_file:
    c = json.load(c_file)

with open('data/community_std_counts/community_std_counts_16s_d.json') as d_file:
    d = json.load(d_file)

with open('./data/rrndb_16s_copies copy.json') as copies_file:
    copy_numbers = json.load(copies_file)

app = dash.Dash(__name__)

app.layout = html.Div([
    html.Div([
        dcc.Markdown(dedent("""
            ### LOD Linearity
        """))
    ], className='row'),
    html.Div([
        html.Div([
            html.P([
                'A459 Concentration'
            ]),
            dcc.Dropdown(
                id='var-select',
                options=[{'label': label, 'value': val} for label, val in
                    [('1E0', 'a'), ('1E2', 'b'), ('1E4', 'c'), ('1E5', 'd')]],
                value='a'
            )
        ], className='two columns'),
        html.Div([
            html.P([
                'Copy Number Adjustment'
            ]),
            daq.BooleanSwitch(
                id='copy-number-switch',
                on=True
            )
        ], className='two columns')
    ], className='row'),
    html.Div([
        dcc.Graph(
            id='graph'
        )
    ], className='row'),
], style={'padding': '10px 50px 50px 50px'})


@app.callback(
    Output('graph', 'figure'),
    [Input('var-select', 'value'),
     Input('copy-number-switch', 'on'),
     ]
)
def update_graph(varSelect, copyNumber):
    if varSelect == 'a':
        count_data = a
    elif varSelect == 'b':
        count_data = b
    elif varSelect == 'c':
        count_data = c
    elif varSelect == 'd':
        count_data = d
    value_dict = filter_data(count_data, copyNumber)
    figure = generate_fig(value_dict)
    return figure


if __name__ == '__main__':
    app.run_server(debug=True, port=9080)
