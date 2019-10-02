import pandas as pd
import os
import plotly.express as px
import numpy as np

def prep_df(df_dict):
    df_prep = pd.DataFrame.from_dict(df_dict)
    df_prep['target_id'] = \
        df['target_id'].apply(lambda x: ('_').join(x.split('_')[-3:-1]))
    df_prep = df_prep.replace({0: np.nan})
    df_prep['sum'] = df_prep.sum(axis=1)
    df_prep['std'] = df_prep.drop(columns='sum').std(axis=1)
    df_prep['mean'] = df_prep.drop(columns=['sum', 'std']).mean(axis=1)
    df_prep['cv'] = df_prep['std'] / df_prep['mean']
    df_prep['samples_observed'] = \
        (~df_prep.drop(columns=['sum', 'std', 'mean', 'cv']).isna()).sum(axis=1)
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
fig_counts = px.scatter(df_counts, x='mean', y='cv',
                        color='samples_observed',
                        hover_data=['mean', 'cv', 'samples_observed', 'target_id'],
                        color_continuous_scale=px.colors.sequential.Rainbow)
df_counts[df_counts['samples_observed'] >= 90].to_csv('kallisto_counts_in_90.csv')
fig_counts.write_html('kallisto_counts.html')
df_tpm = prep_df(tpm_dict)
fig_tpm = px.scatter(df_tpm, x='mean', y='cv',
                     color='samples_observed',
                     hover_data=['mean', 'cv', 'samples_observed', 'target_id'],
                     color_continuous_scale=px.colors.sequential.Rainbow,
                     log_x=True, log_y=True)
df_tpm[df_tpm['samples_observed'] >= 90].to_csv('kallisto_tpm_in_90.csv')
fig_tpm.write_html('kallisto_tpm_log.html')
fig_tpm_lin = px.scatter(df_tpm, x='mean', y='cv',
                         color='samples_observed',
                         hover_data=['mean', 'cv', 'samples_observed', 'target_id'],
                         color_continuous_scale=px.colors.sequential.Rainbow)
fig_tpm_lin.write_html('kallisto_tpm_linear.html')
df_tpm_log_90 = df_tpm.copy()
df_tpm_log_90.loc[(df_tpm['samples_observed'] < 90), ['mean', 'cv']] = np.nan
fig_tpm_90 = px.scatter(df_tpm_log_90, x='mean', y='cv',
                        color='samples_observed',
                        hover_data=['mean', 'cv', 'samples_observed', 'target_id'],
                        color_continuous_scale=px.colors.sequential.Rainbow,
                        log_x=True, log_y=True)
fig_tpm_90.write_html('kallisto_tpm_log_90.html')
