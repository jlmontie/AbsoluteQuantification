import sys
import pandas as pd
import os
import numpy as np
sys.path.insert(0,'../')
import pipeline
import glob
from concurrent.futures import ProcessPoolExecutor, wait
from functools import partial

fqo_path = '../../arup_urine_samples_ge_study/ge_distribution/FastQataloguer_ARUP_Urine_2019-09-11.csv'
fqo = pd.read_csv(fqo_path)
fastq_dir = '/home/jmontgomery/mnt/arup_urine_fqs'
fastq_file_ls = os.listdir(fastq_dir)
pool = ProcessPoolExecutor(20)
###### Missed samples ################
missed = ['190828-1-1-IDBD-D100381-d-06-AHMW75AFXY-GTGCAGACAG-ACGGATGGTA-5d67e5e1_S6_R1_001.postAdapt.postQual.fastq',
'190828-1-1-IDBD-D100382-d-07-AHMW75AFXY-CAATCGGCTG-TTCCTACAGC-5d67e5e1_S7_R1_001.postAdapt.postQual.fastq']
######################################
for fastq_file in missed:
    fastq = os.path.join(fastq_dir, fastq_file)
    accession = '-'.join(fastq.split('-')[3:5])
    outdir = os.path.join('coverage_output', accession)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    index_path_ls = glob.glob("index_*")
    for index_path in index_path_ls:
        output = os.path.join(outdir, os.path.basename(index_path).split('_')[1])
        future = pool.submit(pipeline.main, index_path, fastq, output=output)
