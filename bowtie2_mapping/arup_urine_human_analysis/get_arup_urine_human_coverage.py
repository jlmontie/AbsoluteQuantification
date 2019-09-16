import sys
import pandas as pd
import os
sys.path.insert(0,'../')
import pipeline

index_path = 'index_files'
fqo_path = '../../arup_urine_samples_ge_study/ge_distribution/FastQataloguer_ARUP_Urine_2019-09-11.csv'
fqo = pd.read_csv(fqo_path)

fastq_dir = '/home/jmontgomery/mnt/arup_urine_fqs'
fastq_file_ls = os.listdir(fastq_dir)
pipeline.main(index_path, os.path.join(fastq_dir, fastq_file_ls[2]))
# for fastq in fastq_file_ls:
    # pipeline.main(index_path, os.path.join(fastq_dir, fastq))