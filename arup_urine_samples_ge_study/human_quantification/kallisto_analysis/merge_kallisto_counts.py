import pandas as pd
import os

df = pd.read_csv('kallisto_output/IDBD-D100376/abundance.tsv', sep='\t')
tpm = df['tpm']
kallisto_dir = os.listdir('kallisto_output')
print(kallisto_dir)
kallisto_dir.pop(0)
print(kallisto_dir)
for dir in kallisto_dir:
    df = pd.read_csv(os.path.join('kallisto_output', dir, 'abundance.tsv'),
                     sep='\t')
    print(df.head())
    tpm = tpm + df['tpm']
tpm_combine
print(tpm.head())


