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
    # seq_sple_ls = fqo_merged['Seq Sple']
    accession_ls = fqo_merged['Accession']
    if config['Paths']['SummaryFilePath'] is not None:
        summary_dir_ls = config['Paths']['SummaryFilePath']
        if not isinstance(summary_dir_ls, list):
            summary_dir_ls = [summary_dir_ls] * len(accession_ls)
    else:
        summary_dir_ls = fqo['Diagnostic Output Dir'].tolist()
    dilution_ls = fqo_merged['Dilution Factor'].tolist()
    ctrls_ls = fqo_merged['Control Int Org Names'].tolist()
    org_info_dict = initial_conc.to_dict('index')
    with open(config['rDnaResourceFile']) as rdna_file:
        rDNA_copies = json.load(rdna_file)
    return accession_ls, summary_dir_ls, dilution_ls, ctrls_ls, org_info_dict, rDNA_copies


def fit_model(input_info):
    model = titration_fit(input_info, fit_coverage=config['Fit']['FitCoverage'])
    model.fit()
    print(f"Slope: {model.slope_}")
    print(f"Intercept: {model.intercept_}")
    print(f"Fit metrics: {model.fit_metrics_}")
    model.save_model(outdir=config['Paths']['OutputDir'])
    model.plot_fit(show_fig=config['Output']['ShowPlot'],
                   save_fig=config['Output']['SavePlot'],
                   outdir=config['Paths']['OutputDir'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Train a linear quantification model from a dilution series.')
    parser.add_argument('input_file', help='path to the input yaml file.')
    args = parser.parse_args()
    stream = open(args.input_file, 'r')
    config = yaml.load(stream, Loader=yaml.FullLoader)
    input_info = prep_input_files(config)
    fit_model(input_info)
