import sys
import pandas as pd
# sys.path.insert(0,'../')
# import pipeline
import subprocess

index_path = 'index_files'
fqo_path = '/uufs/chpc.utah.edu/common/home/u0002613/AbsoluteQuantification/AbsoluteQuantification/arup_urine_samples_ge_study/ge_distribution/FastQataloguer_ARUP_Urine_2019-09-11.csv'
fqo = pd.read_csv(fqo_path)
fasta_path = fqo['Original fastq path']
for path in fasta_path:
    subprocess.call(f"cp {path} /scratch/jmontgomery/arup_urine", shell=True)
# fasta_path = fqo['Original fastq path'].values[0]
# pipeline.main(index_path, fasta_path)
