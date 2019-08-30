import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from scipy.signal import savgol_filter as savgol
from scipy import signal
import numpy as np
import os
import re
import glob
import sys
import argparse


def filter_data(coverage_df):
    coverage_df = coverage_df[coverage_df['Coverage'] != 0]
    coverage_df = coverage_df.iloc[::10, :]
    coverage_df = coverage_df.rename(columns={'Coverage': 'Raw'})
    coverage_df_filt = coverage_df.copy()
    q1 = coverage_df_filt['Raw'].quantile(q=0.25)
    q3 = coverage_df_filt['Raw'].quantile(q=0.75)
    iqr = q3 - q1
    upper_threshold = q3 + 4 * iqr
    coverage_df_filt['Clipped'] = coverage_df_filt['Raw'].copy()
    coverage_df_filt.loc[coverage_df_filt['Clipped'] > upper_threshold, 'Clipped'] = np.nan
    b, a = signal.butter(3, 1e-3)
    zi = signal.lfilter_zi(b, a)
    z, _ = signal.lfilter(b, a, coverage_df_filt['Clipped'].dropna(), zi=zi*coverage_df_filt['Clipped'].dropna().values[0])
    z2, _ = signal.lfilter(b, a, z, zi=zi*z[0])
    y = signal.filtfilt(b, a, coverage_df_filt['Clipped'].dropna())
    coverage_df_filt['Filtered'] = coverage_df_filt['Clipped'].copy()
    coverage_df_filt.loc[~coverage_df_filt['Clipped'].isna(), 'Filtered'] = y
    coverage_df_filt_melt = coverage_df_filt.melt(id_vars='Position', value_name='Coverage')
    return y, coverage_df_filt_melt


def generate_fig(coverage_df_filt_melt, coverage_filtered):
    fig = px.line(coverage_df_filt_melt, y='Coverage', color='variable')
    return fig


def process_cov_df(df):
    df = df.rename(columns={0: 'chromosome', 1: 'base', 2: 'Coverage'})
    df['Position'] = df.index.values
    df = df.drop(columns=['chromosome', 'base'])
    return df


def calculate_copy_number(coverage_genomic, coverage_18s):
    df_gen = pd.read_csv(coverage_genomic, sep='\t', header=None)
    df_gen = process_cov_df(df_gen)
    df_18s = pd.read_csv(coverage_18s, sep='\t', header=None)
    df_18s = process_cov_df(df_18s)
    original_len = len(df_gen)
    y_gen, cov_filt_gen = filter_data(df_gen)
    nonzero_len = len(y_gen)
    genome_cov = 100 * nonzero_len / original_len
    y_18s = df_18s.loc[(df_18s['Position'] > 100) & (df_18s['Position'] < 1400), 'Coverage']
    cov_gen_mean, cov_gen_std = np.nanmean(y_gen), np.nanstd(y_gen)
    cov_18s_mean, cov_18s_std = np.nanmean(y_18s), np.nanstd(y_18s)
    copy_num_18s = cov_18s_mean / cov_gen_mean
    copy_num_err = copy_num_18s * np.sqrt((cov_gen_std / cov_gen_mean)**2 + (cov_18s_std / cov_18s_mean)**2)
    str0 = f"Genome coverage: {genome_cov:0.0f}%"
    str1 = f"Mean genome coverage: {cov_gen_mean:0.0f} +/- {cov_gen_std:0.0f}"
    str2 = f"Mean 18S coverage: {cov_18s_mean: 0.0f} +/- {cov_18s_std:0.0f}"
    str3 = f"18S copy number: {copy_num_18s: 0.1f} +/- {copy_num_err:0.1f}"
    print(str0)
    print(str1)
    print(str2)
    print(str3)
    with open(os.path.join(os.getcwd(), 'coverage_stats.txt'), 'w') as stats_out:
        out_str = "value\tmean\tstd\npct_genome_cov\t{:0.0f}\tN/A\ngenome_cov\t{:0.0f}\t{:0.0f}\n18s_cov\t{:0.0f}\t{:0.0f}\ncopy_num\t{:0.0f}\t{:0.0f}"
        stats_out.write(out_str.format(genome_cov, cov_gen_mean, cov_gen_std, cov_18s_mean, cov_18s_std, copy_num_18s, copy_num_err))
    fig = generate_fig(cov_filt_gen, y_gen)
    fig.update_layout(annotations=[dict(x=len(y_gen)/2, y=df_gen['Coverage'].max()*0.75,
                    text=f"{str0}<br>{str1}<br>{str2}<br>{str3}")])
    fig.write_html(os.path.join(os.getcwd(), 'coverage_plot.html'))
    fig.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate 18S copy number from 18S and genomic coverage files.')
    parser.add_argument('-s', '--coverage_18s', help='Path to the 18S coverage file.', required=True)
    parser.add_argument('-g', '--coverage_genomic', help='Path to the genomic coverage file.', required=True)
    args = parser.parse_args()
    calculate_copy_number(args.coverage_genomic, args.coverage_18s)
