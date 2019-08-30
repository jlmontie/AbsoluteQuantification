import os
import mmap
import json
import numpy as np
import argparse
import sys
import pandas as pd
from collections import defaultdict
import re
pd.set_option('mode.chained_assignment', None)

print(os.listdir('mnt'))

def get_counts(groupby_obj):
    read_dict = defaultdict(list)
    conc_dict = defaultdict(list)
    groupby_df = groupby_obj[1]
    summary_path = 'mnt/dna_summries_from_rna_pipe'
    taxid = orgs_taxid[groupby_obj[0]]
    summary_files = os.listdir(summary_path)
    read_count_ls = []
    ctrl_count_ls = []
    conc_ls = []
    seq_sple_series = groupby_df['Seq Sple']
    for seq_sple in seq_sple_series:
        dxsm_bac = [dxsm_path for dxsm_path in summary_files
                    if seq_sple in dxsm_path.lower() and 'bacterial' in dxsm_path]
        if os.stat(os.path.join(summary_path, dxsm_bac[0])).st_size == 0:
            continue
        file = open(os.path.join(summary_path, dxsm_bac[0]))
        s = mmap.mmap(file.fileno(), 0, access=mmap.ACCESS_READ)
        # mmap for checking if taxid is in the file at all.
        if not s.find(f'"taxid": {taxid}'.encode('utf-8')) != -1:
            file.close()
            s.close()
            continue
        dxsm_vir = [dxsm_path for dxsm_path in summary_files
                    if seq_sple in dxsm_path.lower() and 'viral' in dxsm_path]
        conc = groupby_df.loc[groupby_df['Seq Sple'] == seq_sple, 'Dilution'].values[0]
        ctrl_orgs = groupby_df.loc[groupby_df['Seq Sple'] == seq_sple, 'Control Int Org Names'].values[0].split('|')
        alternative_ctrl_orgs = ['enterobacteria phage t7', 'enterobacteria phage pr772']
        for line in file:
            cov_info = json.loads(line)
            if cov_info['taxid'] == taxid:
                gene_info = cov_info['gene_info']
                for gene in gene_info:
                    read_dict[gene['geneid']].append(gene['read_count'])
                    conc_dict[gene['geneid']].append(int(conc))
                    conc_ls.append(int(conc))
        with open(os.path.join(summary_path, dxsm_vir[0])) as file_vir:
            ctrl_counts = []
            for line in file_vir:
                cov_info_vir = json.loads(line)
                if cov_info_vir['name'].lower() in ctrl_orgs:
                    print(ctrl_orgs)
                    print(cov_info_vir['name'].lower())
                    ctrl_counts.append(cov_info_vir['read_count'])
            if len(ctrl_counts) > 0:
                ctrl_count_ls.append(np.mean(np.array(ctrl_counts)).tolist())
        if len(ctrl_counts) < 1:
            print("Searching for alternative ctrls")
            with open(os.path.join(summary_path, dxsm_vir[0])) as file_vir:
                for line in file_vir:
                    cov_info_vir = json.loads(line)
                    if cov_info_vir['name'].lower() in alternative_ctrl_orgs:
                        print("Alternative ctrls found")
                        ctrl_counts.append(cov_info_vir['read_count'])
        if len(ctrl_counts) < 1:
            print("Alternative ctrls not found")
            print(seq_sple)
            print(ctrl_orgs)
        # else:
        #     ctrl_count_ls.append(np.mean(np.array(ctrl_counts)).tolist())
        print(ctrl_counts)
        print("=========================")
        file.close()
    gene_dict = {}
    for gene in read_dict.keys():
        gene_dict.update({
            gene: {
                "Read Counts": read_dict[gene],
                "Concentration": conc_dict[gene]
            }
        })
    return {"Organism": groupby_obj[0], "Gene Counts": gene_dict, "Ctrl Counts": ctrl_count_ls}


qc = pd.read_csv('data/FastQataloguer_FullTable_LOD_DNA_2019-07-17.csv')
orgs_original = [
    'aureus',
    'farcinica',
    'influenzae',
    'pneumoniae',
    'pertussis',
    'wallacei'
]
orgs_replace = {
    's. aureus': 'Saureus',
    'n. farcinica': 'Nfarcinica',
    'bpertussis': 'Bpertussis',
    'hinfluenzae': 'Hinfluenzae',
    'kpneumoniae': 'Kpneumoniae',
    'nwallacei': 'Nwallacei'
}
orgs_taxid = {
    'Hinfluenzae': 727,
    'Bpertussis': 520,
    'Kpneumoniae': 573,
    'Saureus': 1280,
    'Nfarcinica': 37329,
    'Nwallacei': 480035
}

pattern = '|'.join(orgs_original)
qc_orgs = qc[qc['Sample Name'].str.contains(pattern) & ~(qc['Sample Name'].str.contains('hpiv3'))]
qc_orgs['Organism'] = qc_orgs['Sample Name'].str.split('d', expand=True)[0]
print(qc_orgs['Organism'])
qc_orgs['Dilution'] = qc_orgs['Sample Name'].str.split('d', expand=True)[1].str[0].astype(int) * -1 + 1
qc_orgs['Organism'] = qc_orgs['Organism'].str.strip().str.strip()
qc_orgs['Organism'] = qc_orgs['Organism'].replace(orgs_replace)
qc_orgs = qc_orgs.sort_values(['Organism', 'Dilution'])
qc_group = qc_orgs.groupby('Organism')
count_dict_raw = {}
for group in qc_group:
    count_dict_raw.update({f"{orgs_taxid[group[0]]}": get_counts(group)})
with open('data/lod_counts_16s_dna.json', 'w') as outfile:
    json.dump(count_dict_raw, outfile)
