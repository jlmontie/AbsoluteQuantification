import sys
sys.path.insert(0, "../../final_script")
from absoluteQuant import absoluteQuant
import json
import glob
import os
import gzip
import pandas as pd
import numpy as np


def read_summary_files(path):
    summary_object_ls = []
    if path.endswith('.gz'):
        with gzip.open(path,'rt') as summary_file:
            for line in summary_file:
                summary_object_ls.append(json.loads(line))
    else:
        with open(path,'rt') as summary_file:
            for line in summary_file:
                summary_object_ls.append(json.loads(line))
    return summary_object_ls


# summary_file_dir = '/Users/jmontgomery/OneDrive/Documents/IDbyDNA/Code/AbsoluteQuantification/lod_dxsm/idbd_7log_std_titration'
fqo = pd.read_csv('FastQataloguer_Synergy_ZymoII_Titration_2019-09-20.csv')
# fqo['Seq Sple'] = fqo['Seq Sple'].str.upper()
rdna_resource_path = '../../data/rrndb_16s_copies_with_18s.json'
out_dir = 'summary_with_quant'

with open(rdna_resource_path) as resource_file:
    rdna_copy_numbers = json.load(resource_file)

ctrl_ls, seq_sple_ls = [], []
for summary_file_dir, seq_sple in zip(fqo['Diagnostic Output Dir'], fqo['Seq Sple']):
    print(seq_sple)
    # Get summary files
    # summary_file_paths = glob.glob(os.path.join(summary_file_dir, '*' + seq_sple + '*'))
    summary_file_paths = [os.path.join(summary_file_dir, path) for path in os.listdir(summary_file_dir)]
    # print(summary_file_paths)
    summary_file_paths = [path for path in summary_file_paths if seq_sple in path.lower()]
    summary_file_paths = [path for path in summary_file_paths if 'dxsm.out.summary' in path]
    summary_file_paths = [path for path in summary_file_paths if not path.endswith('.done')]
    if len(summary_file_paths) < 3:  # Skip if summary files not present
        continue
    viral_path = [path for path in summary_file_paths if 'rna.viral' in path][0]
    viral_summary = read_summary_files(viral_path)
    bacterial_path = [path for path in summary_file_paths if 'rna.bacterial' in path][0]
    bacterial_summary = read_summary_files(bacterial_path)
    fungpar_path = [path for path in summary_file_paths if 'rna.fungal_parasite' in path][0]
    fungpar_summary = read_summary_files(fungpar_path)
    # Get IC counts
    viral_taxids = [
        10760,
        532076,
        1176767,
        1176765,
        1195074,
        1176434,
        227720,
        1176766,
        1837842,
        482822,
        1527506,
        2053563,
        1075775,
        2079317,
        1075774,
        10759,
        866889,
        1871708
    ]
    viral_counts = []
    for org_info in viral_summary:
        if org_info['taxid'] in viral_taxids:
            viral_counts.append(org_info['read_count'])
    if len(viral_counts) < 1:
        print(f"No IC found for {seq_sple}. Skipping quantification.")
        continue
    # Get quantifications
    bacterial_summary_with_quant = absoluteQuant(viral_counts, bacterial_summary, rdna_copy_numbers, model='synergy') #, quant_mode='read_count')
    fungpar_summary_with_quant = absoluteQuant(viral_counts, fungpar_summary, rdna_copy_numbers, model='synergy') #, quant_mode='read_count')
    # Write modified summary files
    with open(os.path.join(out_dir, os.path.splitext(os.path.basename(bacterial_path))[0]), 'w') as bacterial_out:
        for line in bacterial_summary_with_quant:
            bacterial_out.write(f"{json.dumps(line)}\n")
    with open(os.path.join(out_dir, os.path.splitext(os.path.basename(fungpar_path))[0]), 'w') as fungpar_out:
        for line in fungpar_summary_with_quant:
            fungpar_out.write(f"{json.dumps(line)}\n")
    ctrl_ls.append(np.sum(viral_counts))
    seq_sple_ls.append(seq_sple)
# pd.DataFrame(data={'Seq Sple': seq_sple_ls, 'Ctrl Counts': ctrl_ls}).to_csv('~/Downloads/ctrl_counts.csv')
