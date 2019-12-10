import pandas as pd
import subprocess as sp
import os

fqo = pd.read_csv('../ge_distribution/FastQataloguer_ARUP_Urine_2019-09-11.csv')
fqo_paths = fqo['Original fastq path']
for path in fqo_paths[:5]:
    if not os.path.exists(path):
        print(f"{path} not found. Skipping.")
        continue
    if os.stat(path).st_size == 0:
        print(f"{path} an empty file. Skipping.")
        continue
    sp.call(f"cp {path} ~/mnt/arup_urine_fqs", shell=True)
