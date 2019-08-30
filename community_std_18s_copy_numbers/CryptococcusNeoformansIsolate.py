import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from scipy.signal import savgol_filter as savgol
from scipy import signal
import numpy as np
import os
import re
import glob


def sorted_nicely( l ):
    """ Sort the given iterable in the way that humans expect."""
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)


def generate_fig(coverage_df_filt_melt, coverage_filtered):
    fig = px.line(coverage_df_filt_melt, y='Coverage', color='variable')
    return fig


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


parentdir = 'data/geneious_alignments/cryptococcus_neoformans_full_genome'
cov_file_ls = os.listdir(parentdir)
cov_file_ls = [file for file in cov_file_ls if '18S' not in file]
cov_file_ls = sorted_nicely(cov_file_ls)
filtered_coverage_ls = []
filtered_df_ls = []
chrom_boundaries = []
chrom_start = 0
for cov_file in cov_file_ls:
    print(cov_file)
    cov_df = pd.read_csv(os.path.join(parentdir, cov_file))
    coverage_filtered, filtered_df = filter_data(cov_df)
    chrom_boundaries.append((chrom_start, chrom_start + len(coverage_filtered)))
    chrom_start += len(filtered_df) / 3 # Divide three because 3 columns are melted
    filtered_df_ls.append(filtered_df)
    filtered_coverage_ls.extend(coverage_filtered)
filtered_coverage_mean = np.nanmean(filtered_coverage_ls)
filtered_coverage_std = np.nanstd(filtered_coverage_ls)
print(f"Mean genome coverage: {filtered_coverage_mean:0.0f} +/- {filtered_coverage_std:0.0f}")
fig = generate_fig(pd.concat(filtered_df_ls, ignore_index=True), coverage_filtered)
chrom_shading = []
annotations = []
for idx, (boundaries, cov_file) in enumerate(zip(chrom_boundaries, cov_file_ls)):
    if idx % 2 == 0:
        color = 'rgb(225, 225, 225, 0.5)'
    else:
        color = 'rgb(255, 255, 255, 0)'
    annotation_text = cov_file.split('_')[1].split('.')[0].capitalize()
    # annotation_text.capitalize()
    chrom_shading.append(
        go.layout.Shape(
            type="rect",
            x0=boundaries[0],
            y0=0,
            x1=boundaries[1],
            y1=720,
            fillcolor=color,
            layer='below',
            line=dict(
                color=color,
            )
        )
    )
    annotations.append(
        dict(
            x=np.mean(boundaries),
            y=500,
            text=annotation_text,
            showarrow=False
        )
    )

fig.update_layout(
    shapes=chrom_shading,
    template='none',
    title=f"Mean coverage = {filtered_coverage_mean:0.0f} +/- {filtered_coverage_std:0.0f}",
    xaxis=dict(
        showticklabels=False
    ),
    annotations=annotations
)
fig.show()

cov_18s = pd.read_csv(glob.glob(os.path.join(parentdir, '*_18S.csv'))[0])
cov_18s_mean = cov_18s['Coverage'].mean()
copy_number_18S = cov_18s_mean / filtered_coverage_mean
print(f"Mean 18S coverage: {cov_18s_mean:0.0f}")
print(f"18S copy number: {copy_number_18S:0.1f}")
