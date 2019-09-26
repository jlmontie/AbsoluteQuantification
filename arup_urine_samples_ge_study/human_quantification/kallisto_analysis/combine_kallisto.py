import pandas as pd
import os
import plotly.express as px

def prep_df(df_dict):
    df_prep = pd.DataFrame.from_dict(df_dict)
    df_prep['target_id'] = df['target_id'].apply(lambda x: ('_').join(x.split('_')[-3:-1]))
    df_prep['sum'] = df_prep.sum(axis=1)
    df_prep = df_prep.sort_values('sum', ascending=False)
    df_prep['samples_observed'] = (df_prep.drop(columns='sum') != 0).sum(axis=1)
    df_prep = df_prep.sort_values('samples_observed', ascending=False)
    return df_prep


kallisto_dir = 'kallisto_output'
output_dirs = os.listdir(kallisto_dir)
counts_dict = {}
tpm_dict = {}
i=0
for dir in output_dirs:
    i += 1
    counts = os.path.join(kallisto_dir, dir, 'abundance.tsv')
    df = pd.read_csv(counts, sep='\t')
    counts_dict.update({dir: df['est_counts']})
    tpm_dict.update({dir: df['tpm']})
print(i)

df_counts = prep_df(counts_dict)
fig_counts = px.scatter(df_counts, y='sum',
                        color='samples_observed',
                        hover_data=['sum', 'samples_observed', 'target_id'],
                        color_continuous_scale=px.colors.sequential.Rainbow)
df_counts[df_counts['samples_observed'] > 90].to_csv('kallisto_counts_in_90.csv')
fig_counts.write_html('kallisto_counts.html')
df_tpm = prep_df(tpm_dict)
fig_tpm = px.line(df_counts, y='sum')
fig_tpm.write_html('kallisto_tpm.html')
