import pandas as pd
import yaml
import numpy as np

fqo = pd.read_csv('../../arup_urine_samples_ge_study/ge_distribution/FastQataloguer_ARUP_Urine_2019-10-03_2.csv')
fqo['accession_id'] = fqo['Original fastq path'].apply(lambda x: x.split('/')[-1].split('.')[0].split('_')[0])
random_5 = np.random.choice(np.arange(len(fqo)), 5)
yaml_ls = []
with open('input.yaml', 'w') as yaml_out:
    for idx, row in fqo.iterrows():
        # row = fqo.iloc[num, :]
        fastq_path = row['Original fastq path']
        accession_id = row['accession_id']
        output_dir = '/uufs/chpc.utah.edu/common/home/u0002613/synergy_validation'
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
