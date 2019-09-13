import pandas as pd
import glob
import os
import sys
import json
import gzip
sys.path.insert(0, '../../final_script')
from absoluteQuant import absoluteQuant


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


def get_org_quants(summary_with_quant):
    org_quants = []
    for org_info in summary_with_quant:
        if org_info['taxid'] not in list(organism_taxids.values()):
            continue
        for gene in org_info['gene_info']:
            if gene['geneid'] == 0:
                coverage = gene['coverage']
                if coverage >= 0.97:
                    org_quants.append({
                        'taxid': org_info['taxid'],
                        'name': org_info['name'],
                        'quant': org_info['absolute_quant']
                    })
    return org_quants


taxid_reporting_name_info = {}  # {"taxid"}->{reporting_id:, reporting_name:, compound_id:, class_type:, subclass:}
f = open('explify_reporting_name_info_table.txt')
f.readline()  # skip header
for line in f:
    data = line.strip().split("\t")
    for taxid in data[2].split(","):
        taxid_reporting_name_info[taxid] = {
            "reporting_id": data[1],
            "reporting_name": data[0],
            "compound_id": data[3],
            "class_type": data[4],
            "subclass": [x for x in data[5].split(",")],
            "nucleic_acid": data[6],
            "RNA": data[7],
            "DNA": data[8]
        }

        if taxid_reporting_name_info[taxid]["class_type"] == "viral":  # overload viral dbs tag
            taxid_reporting_name_info[taxid]["DNA"] = taxid_reporting_name_info[taxid]["RNA"]

        if taxid_reporting_name_info[taxid]["RNA"] == ".":
            del taxid_reporting_name_info[taxid]["RNA"]
        if taxid_reporting_name_info[taxid]["DNA"] == ".":
            del taxid_reporting_name_info[taxid]["DNA"]
f.close()
ic_taxids = [
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
with open('organism_taxids.json') as orgs_file:
    organism_taxids = json.load(orgs_file)
rdna_resource_path = '../../data/rrndb_16s_18s.json'
with open(rdna_resource_path) as resource_file:
    rdna_copy_numbers = json.load(resource_file)
out_dir = 'summary_with_quant'
if not os.path.exists(out_dir):
    os.mkdir(out_dir)
summary_dir = '/mnt/ahara_results/synergy_patient_sample_results'
sample_info = pd.read_csv('Nextera_SampleSheets_All_Reproducibility.csv')
sample_ids = sample_info['Sample_ID'].apply(lambda x: '-'.join(x.split('-')[:5]))

for id in sample_ids:
    print(id)
    org_quants = []
    summary_path_base = os.path.join(summary_dir, id)
    summary_paths = glob.glob(summary_path_base + '.rna.*.dxsm.out.summary.gz')
    viral_path = [path for path in summary_paths if 'viral' in path][0]
    viral_summary = read_summary_files(viral_path)
    bacterial_path = [path for path in summary_paths if 'bacterial' in path][0]
    bacterial_summary = read_summary_files(bacterial_path)
    fungpar_path = [path for path in summary_paths if 'fungal_parasite' in path][0]
    fungpar_summary = read_summary_files(fungpar_path)

    # Get IC counts
    ic_counts = []
    for org_info in viral_summary:
        if org_info['taxid'] in ic_taxids:
            ic_counts.append(org_info['read_count'])
    if len(ic_counts) < 1:
        print(f"No IC found for {accession}. Skipping quantification.")
        continue
    
    # Get quantifications
    bacterial_summary_with_quant = absoluteQuant(ic_counts, bacterial_summary,
        rdna_copy_numbers['16s'])
    bacterial_org_quants = get_org_quants(bacterial_summary_with_quant)
    print(bacterial_org_quants)
    fungpar_summary_with_quant = absoluteQuant(ic_counts, fungpar_summary,
        rdna_copy_numbers['18s'])
    fungpar_org_quants = get_org_quants(fungpar_summary_with_quant)
    print(fungpar_org_quants)
    org_quants = bacterial_org_quants + fungpar_org_quants
    print(org_quants)
    with open(os.path.join('simplified_quant', id + '.quant.txt'), 'w') as outfile:
        for line in org_quants:
            outfile.write(f"{json.dumps(line)}\n")

    # Write modified summary files
    with open(os.path.join(out_dir, os.path.splitext(os.path.basename(bacterial_path))[0]), 'w') as bacterial_out:
        for line in bacterial_summary_with_quant:
            bacterial_out.write(f"{json.dumps(line)}\n")
    with open(os.path.join(out_dir, os.path.splitext(os.path.basename(fungpar_path))[0]), 'w') as fungpar_out:
        for line in fungpar_summary_with_quant:
            fungpar_out.write(f"{json.dumps(line)}\n")