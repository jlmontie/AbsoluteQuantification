import os
import json
import plotly
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import numpy as np

summary_parentdir = '/uufs/chpc.utah.edu/common/home/u0002613/AbsoluteQuantification/AbsoluteQuantification/data/arup_urine_summary_files_with_quantification'
outdir = '/uufs/chpc.utah.edu/common/home/u0002613/AbsoluteQuantification/AbsoluteQuantification/arup_urine_samples_ge_study/ge_distribution'
summary_ls = os.listdir(summary_parentdir)
ge_ls = []
accession_ls = []
for summary in summary_ls:
    accession = '-'.join(summary.split('-')[:5])
    with open(os.path.join(summary_parentdir, summary)) as file:
        for line in file:
            cov_info = json.loads(line)
            if cov_info['taxid'] == 562:
                ge_ls.append(cov_info['absolute_quant'])
                accession_ls.append(accession)
accession
ge_log = np.log10(ge_ls)
ge_df = pd.DataFrame(data={'accession': accession_ls, 'log(Genomic Equivalents / ml)': ge_log})
ge_df.to_csv(os.path.join(outdir, 'quantifications.csv'))
fig = px.histogram(ge_df, x='log(Genomic Equivalents / ml)', nbins=10)
fig.update_xaxes(range=[0, 10])
fig.write_html(os.path.join(outdir, 'hist.html'))


