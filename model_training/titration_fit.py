import os
import mmap
import json
import numpy as np
import sys
import pandas as pd
import gzip
from collections import defaultdict
import re
import plotly.express as px
import plotly.graph_objects as go
pd.set_option('mode.chained_assignment', None)
import yaml


class titration_fit(object):
    def __init__(self, args):
        """
        Trains a model for absolute quantification from dilution series of a
        quantified organism mix.
        args is a tuple or list of the following:
            seq_sple_ls (list): the sample Seq Sple identifiers.
            summary_dir_ls (list): the directories containing the output summary
            files for each sample.
            dilution_factor_ls (list): the dilution factor for each sample.
            ctrls_ls (list): the control organism names for each sample. **ONLY
            WORKS FOR PHAGE CTRLS CURRENTLY (i.e., only parses viral summary).
            org_names (dict): names of the organisms in the mix keyed by taxid.
            org_conc (dict): stock concentration of the organisms keyed by
            taxid.
        """
        self.seq_sple_ls = args[0]
        self.summary_dir_ls = args[1]
        self.dilution_ls = args[2]
        self.ctrls_ls = args[3]
        self.org_info = args[4]
        org_counts = []
        print("Calculating organism coverages")
        for seq_sple, summary_dir, dilution, ctrl in \
            zip(self.seq_sple_ls, self.summary_dir_ls, self.dilution_ls,
                self.ctrls_ls):
            org_count_dict_ls = \
                self._get_counts(seq_sple, summary_dir, dilution, ctrl)
            org_counts.extend(org_count_dict_ls)
        print("Combining coverage dictionaries")
        self.org_counts = self._combine_dictionaries(org_counts)


    def _calculate_lower_quart_cov(self, cov_str):
        cov_arr = np.array([int(i) for i in cov_str.split(',')[:-1]])
        # cov_arr = cov_arr[cov_arr > 0]  # Overestimates quant at low coverage
        lower_quart = np.quantile(cov_arr, 0.25)
        lower_quart_cov = cov_arr[cov_arr <= lower_quart]
        return np.mean(lower_quart_cov)


    def _get_ctrl_counts(self, file_vir, ctrl_orgs):
        ctrl_count_ls = []
        for line in file_vir:
            cov_info_vir = json.loads(line)
            if cov_info_vir['name'].lower() in ctrl_orgs:
                ctrl_count_ls.append(cov_info_vir['read_count'])
        if len(ctrl_count_ls) > 0:
            ctrl_count = np.mean(ctrl_count_ls)
        else:
            ctrl_count = np.nan
        return ctrl_count


    def _get_organism_counts(self, file, org_count_dict_ls, dilution, ctrl_count):
        for line in file:
            cov_info = json.loads(line)
            if cov_info['taxid'] in self.org_info:
                taxid = cov_info['taxid']
                organism = self.org_info[taxid]['organism']
                conc = self.org_info[taxid]['stock_concentration'] / dilution
                gene_info = cov_info['gene_info']
                for gene in gene_info:
                    if gene['geneid'] == 0:
                        cov_str = gene['coverage_string']
                        lower_quart_cov = self._calculate_lower_quart_cov(cov_str)
                org_count_dict_ls.extend([{
                    "taxid": taxid,
                    "Organism": organism,
                    "Coverage": lower_quart_cov,
                    "Concentration": conc,
                    "Dilution": dilution,
                    "Ctrl Counts": ctrl_count
                }])
        return org_count_dict_ls


    def _get_summary_file(self, seq_sple, summary_dir, organism):
        summary_files = os.listdir(summary_dir)
        seq_sple_filter = [file for file in summary_files if seq_sple in file]
        organism_filter = [file for file in seq_sple_filter if organism in file]
        lib_filter = [file for file in organism_filter if 'rna' in file]
        summary_filter = [file for file in lib_filter if 'dxsm.out.summary' in file]
        done_filter = [file for file in summary_filter if not file.endswith('done')]
        final_path = os.path.join(summary_dir, done_filter[0])
        # Skip if file is empty
        if os.stat(final_path).st_size == 0:
            return
        if final_path.endswith('.gz'):
            summary_file = gzip.open(final_path,'rt')
        else:
            summary_file = open(final_path,'rt')
        return summary_file


    def _combine_dictionaries(self, dict_ls):
        # Combine list of dictionaries into a single dictionary
        d = defaultdict(lambda: defaultdict(list))
        for count_dict in dict_ls:
            for key in ["Coverage", "Concentration", "Dilution", "Ctrl Counts"]:
                d[count_dict['taxid']][key].append(count_dict[key])
        return d


    def _get_counts(self, seq_sple, summary_dir, dilution, ctrl):
        file_bac = self._get_summary_file(seq_sple, summary_dir, 'bacterial')
        file_fungpar = self._get_summary_file(seq_sple, summary_dir, 'fungal_parasite')
        file_vir = self._get_summary_file(seq_sple, summary_dir, 'viral')

        # Get control counts
        ctrl_orgs = ctrl.split('|')
        ctrl_count = self._get_ctrl_counts(file_vir, ctrl_orgs)
        # Skip if controls not found
        if ctrl_count == np.nan:
            print("Controls not found.")
            print("Sample without controls:")
            print(seq_sple)
            return

        # Get 16S gene counts for bacteria
        org_count_dict_ls = []
        org_count_dict_ls = \
            self._get_organism_counts(file_bac, org_count_dict_ls,
                                      dilution, ctrl_count)

        # Get 18s gene counts for fungal
        org_count_dict_ls = \
            self._get_organism_counts(file_fungpar, org_count_dict_ls,
                                      dilution, ctrl_count)
        file_bac.close()
        file_fungpar.close()
        file_vir.close()
        return org_count_dict_ls
