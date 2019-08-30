import os
import mmap
import json
import numpy as np
import sys
import pandas as pd
from collections import defaultdict
import re
import plotly.express as px
import plotly.graph_objects as go
pd.set_option('mode.chained_assignment', None)


def calculate_lower_quart_cov(cov_str):
    cov_arr = np.array([int(i) for i in cov_str.split(',')[:-1]])
    # cov_arr = cov_arr[cov_arr > 0]  # Overestimates quant at low coverage
    lower_quart = np.quantile(cov_arr, 0.25)
    lower_quart_cov = cov_arr[cov_arr <= lower_quart]
    return np.mean(lower_quart_cov)


def get_ctrl_counts(file_vir, ctrl_orgs):
    ctrl_count_ls = []
    for line in file_vir:
        cov_info_vir = json.loads(line)
        if cov_info_vir['name'].lower() in ctrl_orgs:
            ctrl_count_ls.append(cov_info_vir['read_count'])
    if len(ctrl_count_ls) > 0:
        ctrl_count = np.mean(ctrl_count_ls)
    else:
        ctrl_count = np.nan
    return ctrl_count


def get_organism_counts(file, org_count_dict_ls, dilution, ctrl_count):
    for line in file:
        cov_info = json.loads(line)
        if cov_info['taxid'] in org_taxids:
            taxid = cov_info['taxid']
            organism = panel_orgs[str(taxid)]
            conc = initial_conc[taxid] + dilution
            gene_info = cov_info['gene_info']
            for gene in gene_info:
                if gene['geneid'] == 0:
                    cov_str = gene['coverage_string']
                    lower_quart_cov = calculate_lower_quart_cov(cov_str)
            org_count_dict_ls.extend([{
                "taxid": taxid,
                "Organism": organism,
                "Read Counts": lower_quart_cov,
                "Concentration": conc,
                "Dilution": dilution,
                "Ctrl Counts": ctrl_count
            }])
    return org_count_dict_ls


def get_summary_path(row):
    classification_outdir = row['Diagnostic Output Dir']


def get_counts(row):
    summary_path = 'lod_dxsm/community_std_titration'
    summary_files = os.listdir(summary_path)
    seq_sple = row['Seq Sple']
    dilution = row['Dilution']
    dxsm_bac = [dxsm_path for dxsm_path in summary_files
                if seq_sple in dxsm_path.lower() and 'bacterial' in dxsm_path]
    dxsm_fungpar = [dxsm_path for dxsm_path in summary_files
                if seq_sple in dxsm_path.lower() and 'fungal_parasite' in dxsm_path]
    dxsm_vir = [dxsm_path for dxsm_path in summary_files
                if seq_sple in dxsm_path.lower() and 'viral' in dxsm_path]
    # Skip if file is empty
    if os.stat(os.path.join(summary_path, dxsm_bac[0])).st_size == 0:
        return
    file_bac = open(os.path.join(summary_path, dxsm_bac[0]))
    file_fungpar = open(os.path.join(summary_path, dxsm_fungpar[0]))
    file_vir = open(os.path.join(summary_path, dxsm_vir[0]))

    # Get control counts
    ctrl_orgs = row['Control Int Org Names'].split('|')
    ctrl_count = get_ctrl_counts(file_vir, ctrl_orgs)
    # Skip if controls not found
    if ctrl_count == np.nan:
        print("Controls not found.")
        print("Sample without controls:")
        print(seq_sple)
        return

    # Get 16S gene counts for bacteria
    org_count_dict_ls = []
    org_count_dict_ls = get_organism_counts(file_bac, org_count_dict_ls, dilution, ctrl_count)

    # Get 18s gene counts for fungal
    org_count_dict_ls = get_organism_counts(file_fungpar, org_count_dict_ls, dilution, ctrl_count)
    file_bac.close()
    file_fungpar.close()
    file_vir.close()
    return org_count_dict_ls


fqo = pd.read_csv('model_training/FastQataloguer_CommunityStdTitration_2019-07-26.csv')
concentrations = pd.read_csv('model_training/CommunityStandardConcentrations.csv')
stock_conc = pd.Series(concentrations['stock copies'].values, index=concentrations['taxid']).to_dict()

with open('data/community_standard_panel.json') as panel_file:
    panel_orgs = json.load(panel_file)
org_taxids = [int(key) for key in panel_orgs]

fqo = fqo.sort_values(['Dilution Factor'])
dict_ls = []
for idx, row in fqo.iterrows():
    dict_ls.append(get_counts(row))
# Combine list of dictionaries into a single dictionary
d = defaultdict(lambda: defaultdict(list))
for ls in dict_ls:
    for count_dict in ls:
        for key in ["Read Counts", "Concentration", "Dilution", "Ctrl Counts"]:
            d[count_dict['taxid']][key].append(count_dict[key])
with open(f'model_training/community_std_counts.json', 'w') as outfile:
    json.dump(d, outfile)

