import pandas as pd
from collections import defaultdict
import yaml
from titration_fit import titration_fit
import json
import argparse


def prep_input_files(config):
    fqo = pd.read_csv(config['FastQataloguerOutput'])
    initial_conc = pd.read_csv(config['OrganismStockConcentrations'],
                               index_col='taxid')
    dilution_factors_df = pd.read_csv(config['SampleDilutionFactors'])
    fqo_merged = fqo.merge(dilution_factors_df, on='Accession')
    seq_sple_ls = fqo_merged['Seq Sple']
    if config['Options']['SummaryFilePath'] is not None:
        summary_dir_ls = config['Options']['SummaryFilePath']
        if not isinstance(summary_dir_ls, list):
            summary_dir_ls = [summary_dir_ls] * len(seq_sple_ls)
    else:
        summary_dir_ls = fqo['Diagnostic Output Dir'].tolist()
    dilution_ls = fqo_merged['Dilution Factor'].tolist()
    ctrls_ls = fqo_merged['Control Int Org Names'].tolist()
    org_info_dict = initial_conc.to_dict('index')
    return seq_sple_ls, summary_dir_ls, dilution_ls, ctrls_ls, org_info_dict


def fit_model(input_info):
    model = titration_fit(input_info)


# with open('../data/community_standard_panel.json') as panel_file:
#     panel_orgs = json.load(panel_file)
# org_taxids = [int(key) for key in panel_orgs]

# fqo = fqo.sort_values(['Dilution Factor'])
# dict_ls = []
# for idx, row in fqo.iterrows():
#     dict_ls.append(get_counts(row, concentrations, org_taxids,
#                    summary_path='../lod_dxsm/community_std_titration'))
# # Combine list of dictionaries into a single dictionary
# d = defaultdict(lambda: defaultdict(list))
# for ls in dict_ls:
#     for count_dict in ls:
#         for key in ["Read Counts", "Concentration", "Dilution", "Ctrl Counts"]:
#             d[count_dict['taxid']][key].append(count_dict[key])
# # with open(f'model_training/community_std_counts.json', 'w') as outfile:
# #     json.dump(d, outfile)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Train a linear quantification model from a dilution series.')
    parser.add_argument('input_file', help='path to the input yaml file.')
    args = parser.parse_args()
    stream = open(args.input_file, 'r')
    config = yaml.load(stream)
    input_info = prep_input_files(config)
    fit_model(input_info)
