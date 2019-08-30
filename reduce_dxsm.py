from post_summary_blast import reduce_summary
import os
import pandas as pd
import subprocess as sp


def process_dxsm(dxsm_path, outpath):
    blast_dir = '~/OneDrive/Documents/Code/Python/Bioinformatics/ncbi-blast-2.9.0+/bin/'
    reduce_summary(dxsm_path, blast_dir, outpath)


def bulk_process_dxsm(parentdir, summary_files):
    for idx, file in enumerate(summary_files):
        print(f"Processing file {idx}: {file}")
        if file.endswith('.gz'):
             p1 = sp.Popen(['gunzip', os.path.join(parentdir, file)])
             p1.wait()
             file = file.strip('.gz')
        in_file = os.path.join(parentdir, file)
        output_dir = 'lod_dxsm/reduced_dxsm/'
        out_file = os.path.join(output_dir, file)
        process_dxsm(in_file, out_file)


dxsm_parentdir = 'lod_dxsm/LOD_summary_files/'
lod_files = os.listdir(dxsm_parentdir)
# temp line below to prevent reprocessing bacterial dxsm summaries
bac_vir_summaries = [file for file in lod_files if
                     file.endswith('rna.bacterial.dxsm.out.summary') or
                     file.endswith('rna.bacterial.dxsm.out.summary.gz') or
                     file.endswith('rna.viral.dxsm.out.summary') or
                     file.endswith('rna.viral.dxsm.out.summary.gz')]
bulk_process_dxsm(dxsm_parentdir, bac_vir_summaries)
