import pandas as pd
import os
# from concurrent.futures import ThreadPoolExecutor
# import docker
import subprocess as sp


def kallisto_map(file):
    index = 'GCF_000001405.39_GRCh38.p13_rna_from_genomic.idx'
    outdir = '-'.join(file.split('-')[3:5])
    cmd = f"-v $PWD:/mnt -v {fastq_dir}:/mnt2 kallisto-img kallisto quant --single -l 300 -s 30 -i /mnt/{index} -o /mnt/{outdir} /mnt2/{file}"
    # container = client.containers.run('kallisto-img', command=cmd, detach=True)
    # logs = container.logs()
    # for line in container.logs(stream=True):
    #     print (line.strip())
    sp.call(f"sudo docker run --rm {cmd}", shell=True)


fqo_path = '../../ge_distribution/FastQataloguer_ARUP_Urine_2019-09-11.csv'
fqo = pd.read_csv(fqo_path)
fastq_dir = '/home/jmontgomery/mnt/arup_urine_fqs'
fastq_file_ls = os.listdir(fastq_dir)
# # We want 46 threads in the pool
# pool = ThreadPoolExecutor(46)

# client = docker.from_env()
for file in fastq_file_ls:
    # future = pool.submit(kallisto_map, file, client)
    kallisto_map(file)

