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

summary_file_dir = '/Users/jmontgomery/Desktop/tmp_summary'
sample_info = pd.read_excel('/Users/jmontgomery/OneDrive/Documents/IDbyDNA/Code/AbsoluteQuantification/arup_urine_samples_ge_study/ge_distribution/190904_Urine_Sample_Processing_Log.xlsx')
rdna_resource_path = '/Users/jmontgomery/OneDrive/Documents/IDbyDNA/Code/AbsoluteQuantification/data/rrndb_16s_copies.json'
out_dir = '/Users/jmontgomery/Desktop/tmp_summary_quant_no_t7_expansion'
out_dir_ic_cov = '/Users/jmontgomery/Desktop/tmp_summary_quant_cov_ic'
if not os.path.exists(out_dir_ic_cov):
    os.mkdir(out_dir_ic_cov)

with open(rdna_resource_path) as resource_file:
    rdna_copy_numbers = json.load(resource_file)
# for library in batch['libraries']:
ctrl_ls, accession_ls = [], []
for idx, row in sample_info.iterrows():
    accession = row['Accession #']
    print(accession)
    # Get summary files
    summary_file_paths = glob.glob(os.path.join(summary_file_dir, '*' + accession + '*'))
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
    # # Get viral counts
    viral_taxids = [
        10760
    ]
    viral_counts = []
    for org_info in viral_summary:
        if org_info['taxid'] in viral_taxids:
            viral_counts.append(org_info['read_count'])
    if len(viral_counts) < 1:
        print(f"No IC found for {accession}. Skipping quantification.")
        continue
    # Get quantifications
    bacterial_summary_with_quant, bac_cov_ic_ls = absoluteQuant(viral_counts, bacterial_summary, rdna_copy_numbers)
    fungpar_summary_with_quant, fungpar_cov_ic_ls = absoluteQuant(viral_counts, fungpar_summary, rdna_copy_numbers)
    # Write modified summary files
    with open(os.path.join(out_dir, os.path.splitext(os.path.basename(bacterial_path))[0]), 'w') as bacterial_out:
        for line in bacterial_summary_with_quant:
            bacterial_out.write(f"{json.dumps(line)}\n")
    with open(os.path.join(out_dir, os.path.splitext(os.path.basename(fungpar_path))[0]), 'w') as fungpar_out:
        for line in fungpar_summary_with_quant:
            fungpar_out.write(f"{json.dumps(line)}\n")
    with open(os.path.join(out_dir_ic_cov, os.path.splitext(os.path.basename(bacterial_path))[0]), 'w') as bac_cov_ic_out:
        for line in bac_cov_ic_ls:
            bac_cov_ic_out.write(f"{json.dumps(line)}\n")
    ctrl_ls.append(np.sum(viral_counts))
    accession_ls.append(accession)
pd.DataFrame(data={'Accession': accession_ls, 'Ctrl Counts': ctrl_ls}).to_csv('~/Downloads/ctrl_counts.csv')
