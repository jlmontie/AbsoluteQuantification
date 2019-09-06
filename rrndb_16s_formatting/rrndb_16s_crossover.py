import pandas as pd
import numpy as np
import json
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

with open('/Users/jmontgomery/OneDrive/Documents/IDbyDNA/Code/AbsoluteQuantification/data/rrndb_16s_copies.json') as file:
    rrndb = json.load(file)
uti_orgs = pd.read_csv(
    'rrndb_16s_formatting/synergy_uti_organisms.txt', sep='\t')
uti_taxa = uti_orgs['taxid']
uti_names = uti_orgs['name']
copies_ls = []
for val in rrndb.values():
    copies_ls.append(val['copies'])
mean_copies = np.nanmean(copies_ls)
print(f"Mean 16S copies: {mean_copies:0.1f}")
for taxid, name in zip(uti_taxa, uti_names):
    if str(taxid) not in rrndb:
        genus_taxid = get_genus_taxid(taxid)
        if str(taxid) in rrndb:
            print(f"{taxid} {name} in rrndb at genus level.")
        else:
            print(f"{taxid} {name} not found in rrndb.")
