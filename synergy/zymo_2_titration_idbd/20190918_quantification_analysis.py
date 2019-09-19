import json
import numpy as np
import os
import pandas as pd
import sys
from collections import defaultdict
import plotly
import plotly.graph_objects as go
from plotly.colors import DEFAULT_PLOTLY_COLORS as colors


def combine_dictionaries(dict_ls):
    d = defaultdict(lambda: defaultdict(list))
    for dictionary in dict_ls:
        key_ls = list(dictionary.keys())
        key_ls.remove('taxid')
        for key in key_ls:
            d[dictionary['taxid']][key].append(dictionary[key])
    return d


zymo_orgs = pd.read_csv('ZymoStd2ConcentrationsTitrationIDBD.csv')
titration_key = pd.read_csv('20190918_idbd_zymo_titration_key.csv')
taxids = zymo_orgs['taxid'].tolist()
organisms = zymo_orgs['organism'].tolist()
summary_dir = 'summary_with_quant'
quant_dxsm = os.listdir(summary_dir)

dict_ls = []
for summary in quant_dxsm:
    seq_sple_lower = summary.split('.')[0][:-2].lower()
    dilution_factor = \
        titration_key.loc[titration_key['Seq Sple'].str.contains(seq_sple_lower),
                          'Dilution Factor'].values[0]
    concentrations = zymo_orgs.loc[:, str(dilution_factor)]
    summary_path = os.path.join(summary_dir, summary)
    with open(summary_path) as file:
        for line in file:
            cov_info = json.loads(line)
            if cov_info['taxid'] in taxids:
                for gene in cov_info['gene_info']:
                    if gene['geneid'] == 0:
                        coverage = gene['coverage']
                if coverage < 0.97:
                    continue
                org_idx = taxids.index(cov_info['taxid'])
                if cov_info['absolute_quant'] == 0:
                    calculated_conc = cov_info['absolute_quant']
                else:
                    calculated_conc = np.log10(cov_info['absolute_quant'])
                expected_conc = np.log10(concentrations[org_idx])
                dict_ls.append({'taxid': cov_info['taxid'],
                        'organism': organisms[org_idx],
                        'expected concentration': expected_conc,
                        'calculated concentration': calculated_conc
                })

merged_dict = combine_dictionaries(dict_ls)
merged_df = pd.DataFrame.from_dict(merged_dict, orient='index')

color_dict = dict(zip(taxids, colors)) # for consistent colors for each organism
data = []
for taxid in merged_dict.keys():
    if len(merged_dict[taxid]['calculated concentration']) == 2:
        text = ['', merged_dict[taxid]['organism'][0]]
    else:
        text = merged_dict[taxid]['organism'][0]
    if taxid == 562 or taxid == 1280:
        textposition = 'bottom center'
    else:
        textposition = 'top center'
    if merged_dict[taxid]['calculated concentration'][0] == 0:
        continue
    trace = go.Scatter(
        x=merged_dict[taxid]['expected concentration'],
        y=merged_dict[taxid]['calculated concentration'],
        name=merged_dict[taxid]['organism'][0],
        mode='markers+text',
        marker=dict(
            size=12,
            color=color_dict[taxid]
        ),
        text=text,
        textposition=textposition,
        showlegend=False
    )
    data.append(trace)
data.insert(0, go.Scatter(
    x=[0, 10],
    y=[0, 10],
    showlegend=False,
    mode='lines'
))
layout = go.Layout(
    xaxis=dict(
        range=[0, 10],
        title='Expected Conc. (GE/ml)'
    ),
    yaxis=dict(
        range=[0, 10],
        title='Calculated Conc. (GE/ml)'
    ),
    height=900,
    width=900
)
fig = go.Figure(data, layout)
fig.write_html('quantification_by_coverage_nonzero.html')
# fig.show()