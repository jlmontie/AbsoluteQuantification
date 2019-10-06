import pandas as pd
import numpy as np
import json
from collections import defaultdict
from ncbi_taxonomy_utils import ncbi_taxonomy
ncbi = ncbi_taxonomy()

def get_genus_taxid(taxid):
    idx = 0
    rank = None
    while not (rank == 'genus') or not (idx == 5):
        if len(ncbi.get_path(taxid)) < 2:
            break
        taxid = ncbi.get_path(taxid)[1]
        rank = ncbi.get_rank(taxid)
        idx += 1
    return taxid

with open('rrndb_16s_18s.json') as file:
    rrndb = json.load(file)
uti_orgs = pd.read_csv(
    'synergy_uti_organisms.txt', sep='\t')
uti_taxa_dict = defaultdict(list)
uti_names_dict = defaultdict(list)
for taxid in uti_orgs['taxid']:
    path = ncbi.get_path(taxid)
    if 4751 in path:
        uti_taxa_dict['18s'].append(taxid)
        uti_names_dict['18s'].append(ncbi.get_name(taxid))
    if 2 in path:
        uti_taxa_dict['16s'].append(taxid)
        uti_names_dict['16s'].append(ncbi.get_name(taxid))
print(uti_names_dict)
rdna_out =  open('synergy_uti_rdna_copies.txt', 'w')
rdna_out.write("taxid\torganism\trdna_copies\n")
for org_class in rrndb.keys():
    print(org_class)
    if org_class == '16s':
        class_name = 'bacterial'
        uti_taxa = uti_taxa_dict['16s']
        uti_names = uti_names_dict['16s']
    elif org_class == '18s':
        class_name = 'fungal'
        uti_taxa = uti_taxa_dict['18s']
        uti_names = uti_names_dict['18s']
    rrndb_org_class = rrndb[org_class]
    copies_ls = []
    for val in rrndb_org_class.values():
        copies_ls.append(val['copies'])
    mean_copies = np.nanmean(copies_ls)
    print(f"Mean {org_class} copies: {mean_copies:0.1f}")
    for taxid, name in zip(uti_taxa, uti_names):
        if str(taxid) not in rrndb_org_class:
            genus_taxid = get_genus_taxid(taxid)
            if str(taxid) in rrndb_org_class:
                print(f"{taxid} {name} in rrndb at genus level.")
                rdna_out.write(f"{taxid}\t{name}\tmean of genus used\n")
            else:
                print(f"{taxid} {name} not found in rrndb.")
                rdna_out.write(f"{taxid}\t{name}\tmean of all {class_name} organisms used\n")
        if str(taxid) in rrndb_org_class:
            rdna_out.write(f"{taxid}\t{name}\t{rrndb_org_class[str(taxid)]['copies']}\n")

