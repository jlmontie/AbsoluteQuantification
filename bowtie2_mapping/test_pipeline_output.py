import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from scipy.signal import savgol_filter as savgol
from scipy import signal
import numpy as np
import os
import re
import glob


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


df = pd.read_csv('coverage_out', sep='\t', header=None)
df = df.rename(columns={0: 'chromosome', 1: 'base', 2: 'Coverage'})
df['Position'] = df.index.values
df = df.drop(columns=['chromosome', 'base'])
y, cov_filt = filter_data(df)
filtered_coverage_mean = np.nanmean(y)
filtered_coverage_std = np.nanstd(y)
print(f"Mean genome coverage: {filtered_coverage_mean:0.0f} +/- {filtered_coverage_std:0.0f}")
fig = generate_fig(cov_filt, y)
fig.show()