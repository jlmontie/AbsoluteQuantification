import pandas as pd
import yaml
import numpy as np
import glob
import os

fa_path_ls = glob.glob('../simulated_reads/*.fa')
fa_path_ls = [path for path in fa_path_ls if not path.startswith('5')]
yaml_ls = []
with open('input_in_silico_validation.yaml', 'w') as yaml_out:
    for fa_path in fa_path_ls:
        fastq_path = fa_path
        accession_id = os.path.basename(fa_path).split('.')[0]
        output_dir = '/home/jmontgomery/mnt/synergy_validation'
        yaml_dict = {
            'accession': accession_id,
            'internal_control_reporting_ids': ['26706_10760'],
            'output_dir': output_dir,
            'read_file': fastq_path,
            'tag': 'None',
            'type': 'rna'
        }
        yaml_ls.append(yaml_dict)
    yaml.dump(yaml_ls, yaml_out)
