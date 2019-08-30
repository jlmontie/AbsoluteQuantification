import pandas as pd
import numpy as np
from ncbi_taxonomy_utils import ncbi_taxonomy
import json


def get_genus_species_taxid(taxid):
    path_ls = ncbi.get_path(taxid)
    rank_ls = [ncbi.get_rank(path) for path in path_ls]
    if 'species' in rank_ls:
        species_index = rank_ls.index('species')
        species_taxid = path_ls[species_index]
    else:
        species_taxid = np.nan
    if 'genus' in rank_ls:
        genus_index = rank_ls.index('genus')
        genus_taxid = path_ls[genus_index]
    else:
        genus_taxid = np.nan
    return species_taxid, genus_taxid


ncbi = ncbi_taxonomy()
rrndb = pd.read_csv('data/rrnDB/rrnDB-5.5.tsv', sep='\t')
rrndb[['Species taxid', 'Genus taxid']] = pd.DataFrame(rrndb['NCBI tax id'].apply(get_genus_species_taxid).tolist(), dtype=pd.Int64Dtype())
rrndb['Genus'] = rrndb['Genus taxid'].apply(ncbi.get_name)
rrndb['Species'] = rrndb['Species taxid'].apply(ncbi.get_name)
rrndb_16s_genus = rrndb.groupby(['Genus', 'Genus taxid']).agg({'16S gene count': 'median'}).reset_index()
rrndb_16s_species = rrndb.groupby(['Species', 'Species taxid']).agg({'16S gene count': 'median'}).reset_index()
gene_counts_dict = {}
for idx, row in rrndb_16s_genus.iterrows():
    gene_counts_dict.update({str(row['Genus taxid']): {'copies': round(row['16S gene count'], 1), 'rank': 'genus', 'name': row['Genus']}})
for idx, row in rrndb_16s_species.iterrows():
    gene_counts_dict.update({str(row['Species taxid']): {'copies': round(row['16S gene count'], 1), 'rank': 'genus', 'name': row['Species']}})

with open('data/rrndb_16s_copies.json', 'w') as outfile:
    json.dump(gene_counts_dict, outfile)

rrndb_16s_genus.to_csv('data/dna_16s_copies_genus.csv', index=False)
rrndb_16s_species.to_csv('data/dna_16s_copies_species.csv', index=False)
