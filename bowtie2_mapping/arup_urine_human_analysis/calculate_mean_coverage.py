import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from scipy import signal
from scipy.stats import mstats
import numpy as np
import os
import re
import glob

all_means = open('coverage_stats_all_samples.txt', 'w')
all_means.write('accession\tmean_geomean\tmean_nonzero_geomean\tmedian_geomean\tmedian_nonzero_geomean\n')
cov_dir_ls = os.listdir('coverage_output')
for dir in cov_dir_ls:
    print(dir)
    cov_dir = os.path.join('coverage_output', dir)
    plot_dir = os.path.join(cov_dir, 'coverage_plots')
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    gene_ls, mean_ls, mean_nonzero_ls, median_ls, median_nonzero_ls = [], [], [], [], []
    for cov_path in glob.glob(os.path.join(cov_dir, "*.cov")):
        cov_file = os.path.basename(cov_path).split('.')[0]
        gene_ls.append(cov_file)
        df = pd.read_csv(cov_path, sep='\t', header=None)
        df = df.rename(columns={0: 'accession', 1: 'position', 2: 'coverage'})
        mean = df['coverage'].mean()
        mean_ls.append(mean)
        mean_nonzero = df.loc[(df['coverage'] > 0), 'coverage'].mean()
        mean_nonzero_ls.append(mean_nonzero)
        median = df['coverage'].median()
        median_ls.append(median)
        median_nonzero = df.loc[(df['coverage'] > 0), 'coverage'].median()
        median_nonzero_ls.append(median_nonzero)
        cov_std = df['coverage'].std()
        print(f"Mean genome coverage: {mean:0.0f} +/- {cov_std:0.0f}")
        fig = px.line(df, x='position', y='coverage')
        fig.write_html(f'{os.path.join(plot_dir, cov_file)}.html')

    with open(os.path.join(cov_dir, 'mean_coverage.txt'), 'w') as mean_cov_file:
        mean_cov_file.write("gene\tmean\tmean_nonzero\tmedian\tmedian_nonzero\n")
        for gene, mean, mean_nonzero, median, median_nonzero in \
            zip(gene_ls, mean_ls, mean_nonzero_ls, median_ls, median_nonzero_ls):
            mean_cov_file.write(f"{gene}\t{mean}\t{mean_nonzero}\t{median}\t{median_nonzero}\n")
        mean_geomean = mstats.gmean(np.array(mean_ls))
        mean_nonzero_geomean = mstats.gmean(np.array(mean_nonzero_ls))
        median_geomean = mstats.gmean(np.array(median_ls))
        median_nonzero_geomean = mstats.gmean(np.array(median_nonzero_ls))
        mean_cov_file.write(f"geo_mean\t{mean_geomean}\t{mean_nonzero_geomean}\t{median_geomean}\t{median_nonzero_geomean}")
    all_means.write(f"{dir}\t{mean_geomean}\t{mean_nonzero_geomean}\t{median_geomean}\t{median_nonzero_geomean}\n")
all_means.close()