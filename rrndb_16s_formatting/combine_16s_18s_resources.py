import json

with open('rrndb_16s_copies.json') as rdna_16s_file:
    resource_16s = json.load(rdna_16s_file)

with open('rrndb_18s_copies.json') as rdna_18s_file:
    resource_18s = json.load(rdna_18s_file)

merged_dict = {
    "16s": resource_16s,
    "18s": resource_18s
}

with open('rrndb_16s_18s.json', 'w') as outfile:
    json.dump(merged_dict, outfile)
