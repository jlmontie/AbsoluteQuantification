import pandas as pd
import yaml

fqo = pd.read_csv('../../arup_urine_samples_ge_study/ge_distribution/FastQataloguer_ARUP_Urine_2019-09-11.csv')
yaml_ls = []
for idx, row in fqo.iterrows():
    