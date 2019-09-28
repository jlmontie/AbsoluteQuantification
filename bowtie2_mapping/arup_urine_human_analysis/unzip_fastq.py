import os
import subprocess as sp

fastq_dir = '/home/jmontgomery/mnt/arup_urine_fqs/'
fastq_ls = os.listdir(fastq_dir)
for fastq in fastq_ls:
    if fastq.endswith('.gz'):
        sp.call(f"gunzip {os.path.join(fastq_dir, fastq)}", shell=True)