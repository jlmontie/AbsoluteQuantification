import json

resource_18s = {
    "4932": {"copies": 48, "rank": "species", "name": "Saccharomyces cerevisiae"},
    "5207": {"copies": 48.8, "rank": "species", "name": "Cryptococcus neoformans"}
}
with open('data/rrndb_16s_copies.json') as rdna_16s_file:
    resource_16s = json.load(rdna_16s_file)

merged_dict = {
    "16s": resource_16s,
    "18s": resource_18s
}

with open('data/rrndb_16s_18s.json', 'w') as outfile:
    json.dump(merged_dict, outfile)
