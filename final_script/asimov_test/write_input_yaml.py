import pandas as pd
import yaml
import numpy as np

fqo = pd.read_csv('/home/jmontgomery/AbsoluteQuantification/arup_urine_samples_ge_study/ge_distribution/FastQataloguer_ARUP_Urine_2019-10-03.csv')
fqo['accession_id'] = fqo['Fastq path'].apply(lambda x: x.split('/')[-1].split('.')[0].split('_')[0])
yaml_ls = []
with open('input.yaml', 'w') as yaml_out:
    for idx, row in fqo.iterrows():
        fastq_path = row['Fastq path']
        accession_id = row['accession_id']
        output_dir = '/data/taxonomer2/jmontgomery/synergy_validation/'
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
