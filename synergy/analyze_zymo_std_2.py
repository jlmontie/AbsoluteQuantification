import json
import numpy as np
import os
import pandas as pd
import sys
from collections import defaultdict
import plotly
import plotly.graph_objects as go


def combine_dictionaries(dict_ls):
    d = defaultdict(lambda: defaultdict(list))
    for dictionary in dict_ls:
        key_ls = list(dictionary.keys())
        key_ls.remove('taxid')
        for key in key_ls:
            d[dictionary['taxid']][key].append(dictionary[key])
    return d


zymo_orgs = pd.read_csv('ZymoStd2Concentrations.csv')
taxids = zymo_orgs['taxid'].tolist()
organisms = zymo_orgs['organism'].tolist()
concentrations = zymo_orgs['stock_concentration'].tolist()
quant_dxsm = os.listdir('zymo_2_dxsm_quant')

dict_ls = []
for summary in quant_dxsm:
    summary_path = os.path.join('zymo_2_dxsm_quant', summary)
    with open(summary_path) as file:
        for line in file:
            cov_info = json.loads(line)
            if cov_info['taxid'] in taxids:
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
# print(dict_ls)
merged_dict = combine_dictionaries(dict_ls)
merged_df = pd.DataFrame.from_dict(merged_dict, orient='index')
print(merged_df)

data = []
for taxid in merged_dict.keys():
    if merged_dict[taxid]['calculated concentration'][0] == 0:
        continue
    trace = go.Scatter(
        x=merged_dict[taxid]['expected concentration'],
        y=merged_dict[taxid]['calculated concentration'],
        name=merged_dict[taxid]['organism'][0],
        mode='markers',
        marker=dict(
            size=12
        )
    )
    data.append(trace)
data.append(go.Scatter(
    x=[0, 10],
    y=[0, 10],
    showlegend=False,
    mode='lines'
))
layout = go.Layout(
    xaxis=dict(
        # range=[6, 10],
        title='Expected Conc. (GE/ml)'
    ),
    yaxis=dict(
        # range=[6, 10],
        title='Calculated Conc. (GE/ml)'
    )
)
fig = go.Figure(data, layout)
fig.write_html('quantification_by_coverage_nonzero.html')
fig.show()