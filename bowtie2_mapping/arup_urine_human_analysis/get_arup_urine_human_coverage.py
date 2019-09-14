import sys
import pandas as pd
sys.path.insert(0,'../')
import pipeline

index_path = 'index_files'
fqo = pd.read_csv('../../arup_urine_samples_ge_study/ge_distribution/FastQataloguer_ARUP_Urine_2019-09-11.csv')
fasta_path = fqo['Original fastq path'].values[0]
pipeline.main(index_path, fasta_path)
