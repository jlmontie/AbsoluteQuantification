import json
import numpy as np
import os
import pandas as pd
import sys
from collections import defaultdict
import plotly
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.colors import DEFAULT_PLOTLY_COLORS as colors


def combine_dictionaries(dict_ls, top_level_key):
    d = defaultdict(lambda: defaultdict(list))
    for dictionary in dict_ls:
        key_ls = list(dictionary.keys())
        key_ls.remove(top_level_key)
        for key in key_ls:
            d[dictionary[top_level_key]][key].append(dictionary[key])
    return d


def add_sublots(fig, organism, row=None, col=None, showlegend=True):
    color_dict = dict(zip([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], colors)) # for consistent colors for each dilution
    for dilution in set(merged_dict_organism[organism]['dilution']):
        if dilution > 10:
            continue
        dilution_idx = [i for i, x in enumerate(merged_dict_organism[organism]['dilution']) if x == dilution]
        fig.add_trace(
            go.Scatter(
                x=np.array(merged_dict_organism[organism]['expected concentration'])[dilution_idx],
                y=np.array(merged_dict_organism[organism]['calculated concentration'])[dilution_idx],
                name=str(dilution),
                mode='markers',
                marker=dict(color=color_dict[dilution]),
                legendgroup=str(dilution),
                showlegend=showlegend
            ), row=row, col=col
        )
    fig.update_xaxes(title='Expected Conc. (GE/ml)', range=[3, 10], row=row, col=col)
    fig.update_yaxes(title='Calculated Conc. (GE/ml)', range=[3, 10], row=row, col=col)
    return fig


zymo_orgs = pd.read_csv('ZymoStd2ConcentrationsTitration.csv')
titration_key = pd.read_csv('20190918_synergy_zymo_titration_key.csv')
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
                        'calculated concentration': calculated_conc,
                        'dilution': dilution_factor + 1
                })

merged_dict = combine_dictionaries(dict_ls, 'taxid')
merged_dict_organism = combine_dictionaries(dict_ls, 'organism')
print(merged_dict_organism.keys())

fig = make_subplots(
    rows=2, cols=3,
    subplot_titles=[
        'Pseudomonas aeruginosa',
        'Listeria Monocytogenes',
        'Bacillus subtilis',
        'Saccharomyces cerevisiae',
        'Escherichia coli',
        'Salmonella enterica'
    ]
)
for col in [1, 2, 3]:
    for row in [1, 2]:
        fig.add_trace(
            go.Scatter(
                x=[0, 10],
                y=[0, 10],
                showlegend=False,
                mode='lines',
                line=dict(color='rgb(0.25, 0.25, 0.25)')
            ), row=row, col=col
        )

fig = add_sublots(fig, 'Pseudomonas aeruginosa', row=1, col=1, showlegend=True)
fig = add_sublots(fig, 'Listeria Monocytogenes', row=1, col=2, showlegend=False)
fig = add_sublots(fig, 'Bacillus subtilis', row=1, col=3, showlegend=False)
fig = add_sublots(fig, 'Saccharomyces cerevisiae', row=2, col=1, showlegend=False)
fig = add_sublots(fig, 'Escherichia coli', row=2, col=2, showlegend=False)
fig = add_sublots(fig, 'Salmonella enterica', row=2, col=3, showlegend=False)
fig.write_html('quantification_coverage.html')
