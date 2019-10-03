import pandas as pd
import subprocess as sp

target_sequences = pd.read_csv('selected_targets.txt', sep='\t')
bacterial_fa = '/uufs/chpc.utah.edu/common/home/idbydna-group1/exbox_v2_190404/rna_db_info_100405/16S_sequences_190405.fa'
fungal_fa = '/uufs/chpc.utah.edu/common/home/idbydna-group1/exbox_v2_190404/rna_db_info_100405/18S_sequences.fa'
for idx, row in target_sequences.iterrows():
    seqid = row['seqid']
    if row['class'] == 'bacterial':
        fa = bacterial_fa
    elif row['class'] == 'fungal':
        fa = fungal_fa
    sp.call(f"grep -A1 {seqid} {fa} > ref_seqs/{row['taxid']}_{row['name']}.fa", shell=True)