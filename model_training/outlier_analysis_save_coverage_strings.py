import yaml
import json
import numpy as np
from titration_fit import titration_fit
from train_model import prep_input_files
from ncbi_taxonomy_utils import ncbi_taxonomy


class Coverage(titration_fit):
    # Redefine class method to output coverage strings for the "Coverage" field
    def _get_organism_counts(self, file, org_count_dict_ls, dilution, ctrl_count):
        for line in file:
            cov_info = json.loads(line)
            if cov_info['taxid'] in self._org_info:
                taxid = cov_info['taxid']
                organism = self._org_info[taxid]['organism']
                conc = self._org_info[taxid]['stock_concentration'] / dilution
                gene_info = cov_info['gene_info']
                for gene in gene_info:
                    if gene['geneid'] == 0:
                        cov_str = gene['coverage_string']
                        cov_arr = np.array([int(i) for i in cov_str.split(',')[:-1]]) / ctrl_count
                        cov_arr = cov_arr.tolist()
                        read_count = gene['read_count']
                org_count_dict_ls.extend([{
                    "taxid": taxid,
                    "Organism": organism,
                    "Coverage": cov_arr,
                    "Read Counts": read_count,
                    "Concentration": conc,
                    "Dilution": dilution,
                    "Ctrl Counts": ctrl_count,
                }])
        return org_count_dict_ls


ncbi = ncbi_taxonomy()
stream = open('config_rubi.yaml', 'r')
config = yaml.load(stream, Loader=yaml.FullLoader)
input_info = prep_input_files(config)
model = Coverage(input_info, fit_coverage=config['Fit']['FitCoverage'])
org_cov = model.org_counts
# Save coverage strings for analysis in Dash app on local machine
with open('coverage_strings.json', 'w') as outfile:
    json.dump(org_cov, outfile)
