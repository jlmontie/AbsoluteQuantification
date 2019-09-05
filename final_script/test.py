from absoluteQuant import absoluteQuant
import json
import glob
import os
import gzip


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
batch_json_path = '/Users/jmontgomery/Desktop/tmp_summary/190830B01.json'
rdna_resource_path = '/Users/jmontgomery/OneDrive/Documents/IDbyDNA/Code/AbsoluteQuantification/data/rrndb_16s_copies.json'
out_dir = '/Users/jmontgomery/Desktop/tmp_summary_quant'

# summary_file_dir = '/srv/idbydna-group3/results/idbd_dev/190828_NB551543_0126_AHMW75AFXY/tax'
# batch_json_path = '/srv/idbydna-group3/results/idbd_dev/190828_NB551543_0126_AHMW75AFXY/batch/190827B02.json'
# rdna_resource_path = '/uufs/chpc.utah.edu/common/home/u0002613/AbsoluteQuantification/AbsoluteQuantification/data/rrndb_16s_copies.json'
# out_dir = '/uufs/chpc.utah.edu/common/home/u0002613/AbsoluteQuantification/AbsoluteQuantification/data/arup_urine_summary_files_with_quantification'


with open(batch_json_path) as batch_file:
    batch = json.load(batch_file)
with open(rdna_resource_path) as resource_file:
    rdna_copy_numbers = json.load(resource_file)
for library in batch['libraries']:
    seq_sple = library['seqSple']
    print(seq_sple)
    # Get summary files
    summary_file_paths = glob.glob(os.path.join(summary_file_dir, seq_sple + '*'))
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
    # Get viral counts
    internal_controls = library['internalControls']['organisms'] # if not hardcoded: [ctrl_org['reportingId'] for ctrl_org in internal_controls]
    viral_reportingIds = [
        '26706_10760',
        '26368_532076',
        '42586_1176767',
        '42584_1176765',
        '42575_1195074',
        '42580_1176434',
        '42563_227720',
        '42585_1176766',
        '27076_1837842',
        '26418_482822',
        '27095_1527506',
        '27118_2053563',
        '26377_1075775',
        '27141_2079317',
        '26708_1075774',
        '26700_10759',
        '37655_866889',
        '26353_1871708'
    ]
    viral_counts = []
    for org_info in viral_summary:
        if org_info['reporting_id'] in viral_reportingIds[0]:
            viral_counts.append(org_info['read_count'])
    print(viral_counts)
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
