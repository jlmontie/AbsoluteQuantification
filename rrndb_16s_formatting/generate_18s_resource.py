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
resource_18s = {
    "4932": {"copies": 48, "rank": "species", "name": "Saccharomyces cerevisiae"},
    "5207": {"copies": 48.8, "rank": "species", "name": "Cryptococcus neoformans"}
}
# {taxid: 18s copies}
copies_18s = {
    4932: 48,
    5207: 48.8
}
resource_18s = {}
for taxid in copies_18s.keys():
    species_taxid, genus_taxid = get_genus_species_taxid(taxid)
    species_name, genus_name = ncbi.get_name(species_taxid), ncbi.get_name(genus_taxid)
    resource_18s.update({
        str(species_taxid): {'copies': copies_18s[taxid], 'rank': 'species', 'name': species_name},
        str(genus_taxid): {'copies': copies_18s[taxid], 'rank': 'genus', 'name': genus_name}
    })

with open('rrndb_18s_copies.json', 'w') as outfile:
    json.dump(resource_18s, outfile)

