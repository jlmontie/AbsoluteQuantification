import os
import mmap
import json
import numpy as np
import sys
import pandas as pd
import gzip
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
                "Coverage": lower_quart_cov,
                "Concentration": conc,
                "Dilution": dilution,
                "Ctrl Counts": ctrl_count
            }])
    return org_count_dict_ls


def get_summary_file(row, organism):
    classification_outdir = row['Diagnostic Output Dir']
    seq_sple = row['Seq Sple']
    summary_files = os.listdir(classification_outdir)
    seq_sple_filter = [file for file in summary_files if seq_sple in file]
    organism_filter = [file for file in seq_sple_filter if organism in file]
    lib_filter = [file for file in organism_filter if 'rna' in file]
    summary_filter = [file for file in lib_filter if 'dxsm.out.summary' in file]
    done_filter = [file for file in summary_filter if not file.endswith('done')]
    final_path = os.path.join(classification_outdir, done_filter[0])
    # Skip if file is empty
    if os.stat(final_path).st_size == 0:
        return
    if final_path.endswith('.gz'):
        summary_file = gzip.open(final_path,'rt')
    else:
        summary_file = open(final_path,'rt')
    return summary_file


def get_counts(row, summary_path=None):
    if summary_path is None:
        summary_path = row['Diagnostic Output Dir']
    summary_files = os.listdir(summary_path)
    seq_sple = row['Seq Sple']
    dilution = row['Dilution Factor']
    file_bac = get_summary_file(row, 'bacterial')
    file_fungpar = get_summary_file(row, 'fungal_parasite')
    file_vir = get_summary_file(row, 'viral')

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
