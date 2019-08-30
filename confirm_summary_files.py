import pandas as pd
import os

lod_seq_sple_df = pd.read_csv('./lod_rna_seq_sple.csv')
lod_seq_sple  = lod_seq_sple_df['Seq Sple'].tolist()
lod_summary_files = os.listdir('lod_dxsm/LOD_summary_files')
qc = pd.read_csv('./ExplifySLC_V2pipe_SampleQC.csv', sep='\t', encoding='utf-16')

for seq_sple in lod_seq_sple:
    lod_ls = [file for file in lod_summary_files if seq_sple in file]
    if len(lod_ls) < 1:
        print(f"{seq_sple} not present.")
        continue
    bac_summary_ls = [file for file in lod_ls if
                      file.endswith('rna.bacterial.dxsm.out.summary') or
                      file.endswith('rna.bacterial.dxsm.out.summary.gz')]
    if len(bac_summary_ls) < 1:
        print(f"{seq_sple} bacterial summary not present.")
    vir_summary_ls = [file for file in lod_ls if
                      file.endswith('rna.viral.dxsm.out.summary') or
                      file.endswith('rna.viral.dxsm.out.summary.gz')]
    if len(vir_summary_ls) < 1:
        print(f"{seq_sple} viral summary not present.")

print(qc.loc[qc['Seq Sple'].isin(lod_seq_sple), 'Sple Name'].tolist())
