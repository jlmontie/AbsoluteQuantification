import pandas as pd
import yaml
import numpy as np
import glob
import os

# banks
yaml_outdir = 'synergy_validation_jan_2020'
outdir = '/data/analysis_group1/synergy_validation/quantification/validation_jan_2020/'
simulated_fa_dir = os.path.join(outdir, 'simulated_sequences')
classification_outdir = os.path.join(outdir, 'classification_pipeline_output')
ic_reporting_id = ['27508_10665']

fa_path_ls = glob.glob(os.path.join(simulated_fa_dir, '*.fa'))
fa_path_ls = [path for path in fa_path_ls if not os.path.basename(path).startswith('5')]
yaml_ls = []
with open(os.path.join(yaml_outdir, 'input_in_silico_validation.yaml'), 'w') as yaml_out:
    for fa_path in fa_path_ls:
        fastq_path = fa_path
        accession_id = os.path.basename(fa_path).split('.')[0]
        yaml_dict = {
            'accession': accession_id,
            'internal_control_reporting_ids': ic_reporting_id,
            'output_dir': classification_outdir,
            'read_file': fastq_path,
            'tag': 'None',
            'type': 'rna'
        }
        yaml_ls.append(yaml_dict)
    yaml.dump(yaml_ls, yaml_out)
