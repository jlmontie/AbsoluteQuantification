import pandas as pd
import os

qc = pd.read_csv('data/ExplifySLC_V2pipe_SampleQC.csv', sep='\t', encoding='utf-16')
fqo = pd.read_csv('data/FastQataloguer_FullTable_LOD_2019-06-24.csv')
lod_bool = qc['Sple Name'].str.contains('D[0-9]')
libtype = qc['Lib Type'] == 'DNA'
qc_lod = qc[lod_bool & libtype]
tax_paths_list = qc_lod['Run Fastq'].apply(lambda x: x.split('/')[:8] + ['tax'])
tax_paths_joined = tax_paths_list.apply(lambda x: '/'.join(x))
summary_paths = tax_paths_joined + '/' + qc_lod['Seq Sple'] + '.dna.bacterial.dxsm.out'
summary_paths.name = 'summary_paths'
summary_paths.to_csv('data/lod_dna_summary_paths.csv', index=False, header=True)
seq_sple = qc_lod['Seq Sple'].tolist()
fqo_lod = fqo[fqo['Seq Sple'].isin(seq_sple)]
fqo_lod.to_csv('data/lod_dna_fqo.csv', index=False)