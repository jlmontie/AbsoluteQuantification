import pandas as pd
import os
import json
import mmap
import numpy as np
import plotly
from plotly import graph_objs as go
from scipy.stats import linregress


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
    for idx, org in enumerate(value_dict.keys()):
        org_dict = value_dict[org]
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
                x=-4.5,
                y=np.max([np.max(value_dict[key]['reads_log']) for key in value_dict]) - idx * 0.075,
                xref='x1',
                yref='y1',
                # text=f'R-squared: {org_dict["r_value"]**2:0.2f}<br>Slope: {org_dict["slope"]:0.2f}',
                text=f'Slope: {org_dict["slope"]:0.2f},   R<sup>2</sup>: {org_dict["r_value"]**2:0.2f}',
                showarrow=False,
                font=dict(
                    size=14,
                    color=colors[idx]
                )
            ),
            dict(
                x=-4.5,
                y=np.max([np.max(value_dict[key]['reads_norm_log']) for key in value_dict]) - idx * 0.33,
                xref='x2',
                yref='y2',
                # text=f'R-squared: {org_dict["r_value_norm"]**2:0.2f}<br>Slope: {org_dict["slope_norm"]:0.2f}',
                text=f'Slope: {org_dict["slope_norm"]:0.2f},   R<sup>2</sup>: {org_dict["r_value_norm"]**2:0.2f}',
                showarrow=False,
                font=dict(
                    size=14,
                    color=colors[idx]
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
            title='log(Reads)',
            anchor='x1'
        ),
        yaxis2=dict(
            title='log(IC Normalized Reads)',
            anchor='x2'
        )
    )
    fig = go.Figure(data, layout)
    return fig


def filter_data(count_dict, fit_range, fit_range_norm):
    out_dict = {}
    for org in count_dict.keys():
        org_dict = count_dict[str(org)]
        conc = np.array(org_dict["Concentration"])
        # D1 = 2.9E7 cfu/ml for Bpertussis
        # D1 = 1.84E7 cfu/ml for Kpneumoniae
        # D1 = 4.7E6 cfu/ml for Saureus
        # D1 = 8.0E4 cfu/ml for Nfarcinica
        reads, ctrls = org_dict['Read Counts'], org_dict['Ctrl Counts']
        copy_number_normalizer = copy_numbers[org]['copies']
        print(copy_number_normalizer)
        reads = np.array(reads) / copy_number_normalizer
        reads = np.array(reads)
        reads_nonzero_idx = np.argwhere(reads > 0).reshape(-1)
        reads_nonzero = reads[reads_nonzero_idx]
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

with open('./data/rrndb_16s_copies.json') as copies_file:
    copy_numbers = json.load(copies_file)


def update_graph(varSelect, fit_range, fit_range_norm):
    if varSelect == 'a':
        count_data = a
    elif varSelect == 'b':
        count_data = b
    elif varSelect == 'c':
        count_data = c
    elif varSelect == 'd':
        count_data = d
    value_dict = filter_data(count_data, fit_range, fit_range_norm)
    figure = generate_fig(value_dict)
    return figure


fig = plotly.subplots.make_subplots
