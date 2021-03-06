import os
import json
import plotly
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from plotly.colors import DEFAULT_PLOTLY_COLORS as colors
import pandas as pd
import numpy as np
from scipy.stats import linregress
import re


def search_string(string):
    pattern = '(>)?\d{1,3}(,\d{3})\W\S+\W\S+\W\S+'
    match_ls = []
    start_pos = 0
    for meaningless_idx in range(5):
        match_obj = re.search(pattern, string[start_pos:])
        if match_obj is None:
            break
        match = match_obj.group()
        match_ls.append(match)
        match_end = match_obj.span()[1]
        start_pos = match_end
    match_ls = list(set(match_ls))
    return match_ls


def filter_re_match(match):
    organism_ls = []
    concentration_ls = []
    for index, value in match.iteritems():
        for item in value:
            item_components = item.split(' ')
            concentration = item_components[0]
            concentration = concentration.replace('>', '')
            concentration = concentration.replace(',', '')
            organism = ' '.join(item_components[2:])
            organism_ls.append(organism)
            concentration_ls.append(concentration)
    return organism_ls, concentration_ls


summary_parentdir = '/Users/jmontgomery/Desktop/tmp_summary_quant'
sample_info = pd.read_excel('/Users/jmontgomery/OneDrive/Documents/IDbyDNA/Code/AbsoluteQuantification/arup_urine_samples_ge_study/ge_distribution/190904_Urine_Sample_Processing_Log.xlsx')
outdir = '/Users/jmontgomery/Desktop/'
re_match = sample_info['RESULT LONG TEXT'].apply(search_string)

org_taxids = {
    'Escherichia coli': [562],
    'Beta Hemolytic': [1311],
    'Coagulase negative': [
        1294,
        569857,
        522262,
        1282,
        29388,
        29380,
        33028,
        1292,
        45972,
        1283,
        586733,
        1290,
        1276936,
        28035,
        29379,
        1286,
        1281,
        70255,
        70258,
        170573,
        555791,
        29385,
        246432,
        1293,
        61015,
        1288,
        29382,
        214473,
        29378,
        29384,
        1296,
        150056,
        42858,
        643214,
        71237
    ],
    'Enterococcus faecalis': [1351],
    'Streptococcus dysgalactiae': [1334, 119602],
    'Aerococcus urinae': [1376],
    'Enterococcus species': [1350, 1351, 1352],
    'Klebsiella pneumoniae': [573]
}

sample_info['ORGANISM'], sample_info['CONCENTRATION'] = filter_re_match(re_match)
sample_info['TAXID'] = sample_info['ORGANISM'].copy()
for organism in sample_info['ORGANISM'].unique():
    org_match_idx = sample_info.index[sample_info['ORGANISM'] == organism]
    for idx in org_match_idx:
        sample_info.at[sample_info.index[idx], 'TAXID'] = org_taxids[organism]

summary_ls = os.listdir(summary_parentdir)
ge_ls = []
accession_ls = []
arup_conc_ls = []
plate_conc_ls = []
organism_ls = []
for idx, row in sample_info[~sample_info['Accession #'].isna()].iterrows():
    accession = row['Accession #']
    arup_conc = int(row['CONCENTRATION'])
    org = row['ORGANISM']
    bac_summary = [file for file in summary_ls if accession in file and 'bacterial' in file]
    if len(bac_summary) == 0: # Skip if summary for sample not present
        print(f"{accession} skipped. No IC present.")
        continue
    bac_summary = bac_summary[0]
    with open(os.path.join(summary_parentdir, bac_summary)) as file:
        ge_log = np.nan
        for line in file:
            cov_info = json.loads(line)
            if cov_info['taxid'] in row['TAXID']:
                if cov_info['absolute_quant'] != np.nan and cov_info['absolute_quant'] != 0:
                    ge_log = np.log10(cov_info['absolute_quant'])
                else:
                    print(f"{accession} skipped. Organism not present.")
                    ge_log = np.nan
    accession_ls.append(accession)
    arup_conc_ls.append(np.log10(arup_conc))
    organism_ls.append(org)
    if 'Count' in row:
        plate_conc_ls.append(row['Count'])
    ge_ls.append(ge_log)

plate_conc_arr = np.array(plate_conc_ls)
plate_conc_arr[~np.isnan(plate_conc_arr)] = \
    np.log10(plate_conc_arr[~np.isnan(plate_conc_arr)])
ge_df = pd.DataFrame(data={
    'Accession': accession_ls,
    'Detected Organism': organism_ls,
    'log(Genomic Equivalents / ml)': ge_ls,
    'log(cfu/ml)': arup_conc_ls,
    'Plate count log(cfu/ml)': plate_conc_arr
})
# Manually add quantification for disagreeing detections
# ge_df['Detections'] = 'Concordant'
# ge_df.loc[ge_df['Accession'] == 'IDBD-D100414', 'log(Genomic Equivalents / ml)'] = 7.638295708279279
# ge_df.loc[ge_df['Accession'] == 'IDBD-D100415', 'log(Genomic Equivalents / ml)'] = 7.743110586
# ge_df.loc[ge_df['Accession'] == 'IDBD-D100416', 'log(Genomic Equivalents / ml)'] = 8.434176604
# ge_df.loc[ge_df['Accession'] == 'IDBD-D100414', 'Detections'] = 'Discordant - P.aeruginosa'
# ge_df.loc[ge_df['Accession'] == 'IDBD-D100415', 'Detections'] = 'Discordant - P.aeruginosa'
# ge_df.loc[ge_df['Accession'] == 'IDBD-D100416', 'Detections'] = 'Discordant - P.aeruginosa'

ge_df.to_csv(os.path.join(outdir, 'quantifications.csv'))
ge_df = ge_df[~ge_df['log(Genomic Equivalents / ml)'].isna()]
fig = px.histogram(ge_df, x='log(Genomic Equivalents / ml)', nbins=20)
fig.update_xaxes(range=[4, 9])
fig.write_html(os.path.join(outdir, 'hist.html'))
fig_cor = px.scatter(ge_df, x='log(cfu/ml)', y='log(Genomic Equivalents / ml)')
fig_cor.update_xaxes(range=[4, 5.5], dtick=0.5)
fig_cor.update_yaxes(range=[5, 9])
fig_cor.write_html(os.path.join(outdir, 'corr.html'))

lin_fit_x = ge_df.loc[~(ge_df['Plate count log(cfu/ml)'].isna()),
                      'Plate count log(cfu/ml)']
lin_fit_y = ge_df.loc[~(ge_df['Plate count log(cfu/ml)'].isna()),
                      'log(Genomic Equivalents / ml)']
slope, intercept, rval, pval, stderr = linregress(lin_fit_x, lin_fit_y)
fig_sub = make_subplots(rows=1, cols=2)
for idx, organism in enumerate(ge_df['Detected Organism'].unique()):
    fig_sub.add_trace(
        go.Histogram(
            x=ge_df.loc[ge_df['Detected Organism'] == organism, 'log(Genomic Equivalents / ml)'],
            # nbinsx=20,
            legendgroup=organism,
            name=organism,
            marker=dict(
                color=colors[idx]
            )
        ),
        row=1, col=1
    )
    fig_sub.add_trace(
        go.Scatter(
            x=ge_df.loc[(ge_df['Detected Organism'] == organism) &
                        ~(ge_df['Plate count log(cfu/ml)'].isna()), 'Plate count log(cfu/ml)'],
            y=ge_df.loc[(ge_df['Detected Organism'] == organism) &
                        ~(ge_df['Plate count log(cfu/ml)'].isna()), 'log(Genomic Equivalents / ml)'],
            text=ge_df.loc[(ge_df['Detected Organism'] == organism) &
                        ~(ge_df['Plate count log(cfu/ml)'].isna()), 'Accession'],
            mode='markers',
            showlegend=False,
            name=organism,
            marker=dict(
                color=colors[idx],
                size=12
            ),
            hoverinfo='text'
        ),
        row=1, col=2
    )
    fig_sub.add_trace(
        go.Scatter(
            x=np.linspace(lin_fit_x.min(), lin_fit_x.max()),
            y=slope * np.linspace(lin_fit_x.min(), lin_fit_x.max()) + intercept,
            showlegend=False
        ),
        row=1, col=2
    )
fig_sub.update_xaxes(
        title_text='log(Genomic Equivalents / ml)',
    row=1, col=1
)
fig_sub.update_xaxes(
    title_text='Plate count log(cfu/ml)',
    row=1, col=2
)
fig_sub.update_yaxes(
    title_text='Counts',
    row=1, col=1
)
fig_sub.update_yaxes(
    title_text='log(Genomic Equivalents / ml)',
    row=1, col=2
)
fig_sub.update_layout(annotations=[dict(
    x=4.75,
    y=8,
    text=f"Pearson R:<br>{rval:0.2f}",
    showarrow=False,
    xref='x2',
    yref='y2',
    font=dict(
        size=18
    )
)])
fig_sub.write_html(os.path.join(outdir, 'hist_corr.html'))