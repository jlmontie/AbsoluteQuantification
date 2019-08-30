import pandas as pd
from copy_number_calculator import get_counts

fqo = pd.read_csv('model_training/FastQataloguer_CommunityStdTitration_2019-07-26.csv')
concentrations = pd.read_csv('model_training/CommunityStandardConcentrations.csv')
stock_conc = pd.Series(concentrations['stock copies'].values, index=concentrations['taxid']).to_dict()

with open('data/community_standard_panel.json') as panel_file:
    panel_orgs = json.load(panel_file)
org_taxids = [int(key) for key in panel_orgs]

fqo = fqo.sort_values(['Dilution Factor'])
dict_ls = []
for idx, row in fqo.iterrows():
    dict_ls.append(get_counts(row))
# Combine list of dictionaries into a single dictionary
d = defaultdict(lambda: defaultdict(list))
for ls in dict_ls:
    for count_dict in ls:
        for key in ["Read Counts", "Concentration", "Dilution", "Ctrl Counts"]:
            d[count_dict['taxid']][key].append(count_dict[key])
# with open(f'model_training/community_std_counts.json', 'w') as outfile:
#     json.dump(d, outfile)