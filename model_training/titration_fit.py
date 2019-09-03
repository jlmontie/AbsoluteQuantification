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
from plotly.subplots import make_subplots
import plotly.graph_objects as go
pd.set_option('mode.chained_assignment', None)
import yaml
from scipy.stats import linregress
import numpy as np
from ncbi_taxonomy_utils import ncbi_taxonomy
ncbi = ncbi_taxonomy()


class titration_fit(object):
    def __init__(self, args, fit_coverage=True):
        """
        Trains a model for absolute quantification from dilution series of a
        quantified organism mix.
        args:
            args (tuple or list) containing:
                seq_sple_ls (list): the sample Seq Sple identifiers.
                summary_dir_ls (list): the directories containing the output summary
                files for each sample.
                dilution_factor_ls (list): the dilution factor for each sample.
                ctrls_ls (list): the control organism names for each sample. **ONLY
                WORKS FOR PHAGE CTRLS CURRENTLY (i.e., only parses viral summary).
                org_names (dict): names of the organisms in the mix keyed by taxid.
                org_conc (dict): stock concentration of the organisms keyed by
                taxid.
            fit_coverage (bool): fit model with lower quartile of coverage,
            otherwise fit to rDNA read counts. Default True.
        """
        self._seq_sple_ls = args[0]
        self._summary_dir_ls = args[1]
        self._dilution_ls = args[2]
        self._ctrls_ls = args[3]
        self._org_info = args[4]
        self._rdna_copies = args[5]
        org_counts = []
        print("Calculating organism coverages")
        for seq_sple, summary_dir, dilution, ctrl in \
            zip(self._seq_sple_ls, self._summary_dir_ls, self._dilution_ls,
                self._ctrls_ls):

            org_count_dict_ls = \
                self._get_counts(seq_sple, summary_dir, dilution, ctrl)
            for org_count_dict in org_count_dict_ls:
                org_count_dict.update({'Accession': seq_sple})
            org_counts.extend(org_count_dict_ls)

        print("Combining coverage dictionaries")
        self.org_counts = self._combine_dictionaries(org_counts)
        self._fit_coverage = fit_coverage
        self.slope_ = np.nan
        self.intercept_ = np.nan
        self.fit_metrics_ = {}
        self._coverage_log = []
        self._concentration_log = []


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
            if cov_info['taxid'] in self._org_info:
                taxid = cov_info['taxid']
                organism = self._org_info[taxid]['organism']
                conc = self._org_info[taxid]['stock_concentration'] / dilution
                gene_info = cov_info['gene_info']
                for gene in gene_info:
                    if gene['geneid'] == 0:
                        cov_str = gene['coverage_string']
                        lower_quart_cov = self._calculate_lower_quart_cov(cov_str)
                        read_count = gene['read_count']
                org_count_dict_ls.extend([{
                    "taxid": taxid,
                    "Organism": organism,
                    "Coverage": lower_quart_cov,
                    "Read Counts": read_count,
                    "Concentration": conc,
                    "Dilution": dilution,
                    "Ctrl Counts": ctrl_count,
                }])
        return org_count_dict_ls


    def _get_summary_file(self, seq_sple, summary_dir, organism):
        summary_files = os.listdir(summary_dir)
        seq_sple_filter = [file for file in summary_files if seq_sple in file.lower()]
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
            for key in ["Coverage", "Concentration", "Dilution", "Ctrl Counts",
                        "Read Counts", "Accession"]:
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


    def fit(self):
        """
        Fit the model with the data given in the initialization of the class.
        args:
            None
        """
        if self._fit_coverage:
            cov_key = 'Coverage'
        else:
            cov_key = 'Read Counts'
        cov_norm_log_ls = []
        conc_log_ls = []
        accession_ls = []
        taxids_ls = []
        org_names_ls = []
        for org in self.org_counts.keys():
            org_dict = self.org_counts[org]
            conc = np.array(org_dict["Concentration"])
            reads, ctrls = org_dict[cov_key], org_dict['Ctrl Counts']
            copy_number_normalizer = self._rdna_copies[str(org)]['copies']
            reads = np.array(reads) / copy_number_normalizer
            reads_nonzero_idx = np.argwhere(reads > 0).reshape(-1)
            reads_nonzero = reads[reads_nonzero_idx]
            ctrls_nonzero = np.array(ctrls)[reads_nonzero_idx]
            conc_nonzero = np.log10(conc[reads_nonzero_idx])
            reads_log = np.log10(reads_nonzero)
            slope, intercept, r_value, p_value, std_err = linregress(reads_log, conc_nonzero)
            fit_vals_y = slope * reads_log + intercept
            fit_vals_x = reads_log
            # Normalized
            reads_norm = reads_nonzero / ctrls_nonzero
            reads_norm_log = np.log10(reads_norm)
            slope_norm, intercept_norm, r_value_norm, p_value_norm, std_err_norm = \
                linregress(reads_norm_log, conc_nonzero)
            fit_vals_norm_y = slope_norm * reads_norm_log + intercept_norm
            fit_vals_norm_x = reads_norm_log
            cov_norm_log_ls.extend(reads_norm_log)
            conc_log_ls.extend(conc_nonzero)
            accession_ls.extend(org_dict['Accession'])
            org_names_ls.extend([ncbi.get_name(org)] * len(org_dict['Accession']))
            taxids_ls.extend([org] * len(org_dict['Accession']))
        slope, intercept, rval, pval, stderr = linregress(cov_norm_log_ls, conc_log_ls)
        self.slope_ = slope
        self.intercept_ = intercept
        self.fit_metrics_ = {'R-sqr': rval**2, 'P-value': pval, 'Stderr': stderr}
        self._coverage_log = np.array(cov_norm_log_ls)
        self._concentration_log = np.array(conc_log_ls)
        self._accessions_extended = accession_ls
        self._org_names_extended = org_names_ls
        self._taxids_extended = taxids_ls


    def save_model(self, outdir=None):
        """
        Save the model coefficients.
        args:
            outdir (str): Output directory to save coefficients. File will be
            named 'model_coefficients.txt'.
        """
        if outdir is None:
            outdir = os.getcwd()
        with open(os.path.join(outdir, 'model_coefficients.txt'), 'w') as outfile:
            outfile.write(f"slope\t{self.slope_}\nintercept\t{self.intercept_}")


    def plot_fit(self, outdir=None, show_fig=True, save_fig=True):
        """
        Creat plots for the trained model. Includes two subplots 1) regression
        2) residuals.
        args:
            outdir (str): Output directory to save plot. Only valid with
            save_fig=True.
            show_fig (bool): Display the figure in a browser. Default True.
            save_fig (bool): Save the figure in html format. Default True.
        """
        if outdir is None:
            outdir = os.getcwd()
        if self._fit_coverage:
            xtitle = 'log(Normalized Coverage)'
        else:
            xtitle = 'log(Normalized Read Count)'
        text = [f"Accession: {accession}<br>Organism: {org}<br>Taxid: {taxid}"
                for accession, org, taxid in zip(self._accessions_extended,
                                                 self._org_names_extended,
                                                 self._taxids_extended)]
        y_hat = self.slope_ * self._coverage_log + self.intercept_
        fig = make_subplots(rows=1, cols=2, subplot_titles=('Fit', 'Residuals'))
        # Data scatter
        fig.add_trace(
            go.Scatter(
                x=self._coverage_log,
                y=self._concentration_log,
                mode='markers',
                showlegend=False,
                hovertext=text
            ),
            row=1, col=1
        )
        # Fit line
        fig.add_trace(
            go.Scatter(
                x=self._coverage_log,
                y=y_hat,
                showlegend=False,
                mode='lines'
            ),
            row=1, col=1
        )
        # Residuals
        fig.add_trace(
            go.Histogram(
                x=y_hat - self._concentration_log,
                showlegend=False
            ),
            row=1, col=2
        )
        fig.update_xaxes(title_text=xtitle, row=1, col=1)
        fig.update_yaxes(title_text='log(Genomic Equivalents) / ml', row=1, col=1)
        fig.update_xaxes(title_text='Absolute Deviation', row=1, col=2)
        annotations=[
            dict(
                xref='x1',
                yref='y1',
                x=np.mean(self._coverage_log),
                y=np.max(self._concentration_log),
                text=f"Slope: {self.slope_:0.2f}<br>R-sq: {self.fit_metrics_['R-sqr']:0.2f}<br>Std err: {self.fit_metrics_['Stderr']:0.2f}"
            )
        ]
        if show_fig:
            fig.show()
        if save_fig:
            fig.write_html(os.path.join(outdir, 'regression_plot.html'))
