import pandas as pd
import yaml
import numpy as np

fqo = pd.read_csv('../../arup_urine_samples_ge_study/ge_distribution/FastQataloguer_ARUP_Urine_2019-09-11.csv')
random_5 = np.random.choice(np.arange(len(fqo)), 5)
yaml_ls = []
with open('input.yaml', 'w') as yaml_out:
    for num in random_5:
        row = fqo.iloc[num, :]
        fastq_path = row['Original fastq path']
        seq_sple = '-'.join(fastq_path.split('/')[-1].split('-')[:11])
        output_dir = '/uufs/chpc.utah.edu/common/home/u0002613/AbsoluteQuantification/AbsoluteQuantification/final_script/rubi_test_files/output_files'
        yaml_dict = {
            'accession': seq_sple,
            'internal_control_reporting_ids': "['26706_10760']",
            'output_dir': output_dir,
            'read_file': fastq_path,
            'tag': 'None',
            'type': 'rna'
        }
        yaml_ls.append(yaml_dict)
    yaml.dump(yaml_ls, yaml_out)
