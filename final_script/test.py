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
# batch_json_path = '/Users/jmontgomery/Desktop/tmp_summary/190903B02.json'
sample_info = pd.read_excel('/Users/jmontgomery/OneDrive/Documents/IDbyDNA/Code/AbsoluteQuantification/arup_urine_samples_ge_study/ge_distribution/190904_Urine_Sample_Processing_Log.xlsx')
rdna_resource_path = '/Users/jmontgomery/OneDrive/Documents/IDbyDNA/Code/AbsoluteQuantification/data/rrndb_16s_copies.json'
out_dir = '/Users/jmontgomery/Desktop/tmp_summary_quant'

# summary_file_dir = '/srv/idbydna-group3/results/idbd_dev/190828_NB551543_0126_AHMW75AFXY/tax'
# batch_json_path = '/srv/idbydna-group3/results/idbd_dev/190828_NB551543_0126_AHMW75AFXY/batch/190827B02.json'
# rdna_resource_path = '/uufs/chpc.utah.edu/common/home/u0002613/AbsoluteQuantification/AbsoluteQuantification/data/rrndb_16s_copies.json'
# out_dir = '/uufs/chpc.utah.edu/common/home/u0002613/AbsoluteQuantification/AbsoluteQuantification/data/arup_urine_summary_files_with_quantification'


# with open(batch_json_path) as batch_file:
#     batch = json.load(batch_file)
with open(rdna_resource_path) as resource_file:
    rdna_copy_numbers = json.load(resource_file)
# for library in batch['libraries']:
ctrl_ls, accession_ls = [], []
for idx, row in sample_info.iterrows():
    accession = row['Accession #']
    # seq_sple = library['seqSple']
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
    # internal_controls = library['internalControls']['organisms'] # if not hardcoded: [ctrl_org['reportingId'] for ctrl_org in internal_controls]
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
        print(f"No IC found for {accession}. Skipping quantification.")
        continue
    # Get quantifications
    bacterial_summary_with_quant = absoluteQuant(viral_counts, bacterial_summary, rdna_copy_numbers)
    fungpar_summary_with_quant = absoluteQuant(viral_counts, fungpar_summary, rdna_copy_numbers)
    # Write modified summary files
    with open(os.path.join(out_dir, os.path.splitext(os.path.basename(bacterial_path))[0]), 'w') as bacterial_out:
        for line in bacterial_summary_with_quant:
            bacterial_out.write(f"{json.dumps(line)}\n")
    with open(os.path.join(out_dir, os.path.splitext(os.path.basename(fungpar_path))[0]), 'w') as fungpar_out:
        for line in fungpar_summary_with_quant:
            fungpar_out.write(f"{json.dumps(line)}\n")
    ctrl_ls.append(np.sum(viral_counts))
    accession_ls.append(accession)
pd.DataFrame(data={'Accession': accession_ls, 'Ctrl Counts': ctrl_ls}).to_csv('~/Downloads/ctrl_counts.csv')
