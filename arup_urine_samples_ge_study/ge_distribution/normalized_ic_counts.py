import pandas as pd
import json
import glob
import os
import gzip
import pandas as pd
import numpy as np
import plotly as px


def read_summary_files(path):
    summary_object_ls = []
    if path.endswith('.gz'):
        with gzip.open(path,'rt') as summary_file:
            for line in summary_file:
                summary_object_ls.append(json.loads(line))
    else:
        with open(path,'rt') as summary_file:
            for line in summary_file:
                summary_object_ls.append(json.loads(line))
    return summary_object_ls


def get_read_counts(summaryObjectList, taxid_ls):
    read_count, mean_coverage, median_coverage, q1 = np.nan, np.nan, np.nan, np.nan
    for idx, coverageInfo in enumerate(summaryObjectList):
        taxid = coverageInfo['taxid']
        if taxid in taxid_ls:
            print('organism found')
            for gene in coverageInfo['gene_info']:
                if gene['geneid'] == 0:
                    read_count = gene['read_count']
                    coverageStr = gene['coverage_string']
                    coverage = coverageStr.split(',')[:-1]
                    coverage = np.array([int(item) for item in coverage])
                    mean_coverage = np.mean(coverage)
                    median_coverage = np.median(coverage)
                    q1 = np.quantile(coverage, 0.25)
                    coverageQ1 = coverage[coverage <= q1]
    return read_count, mean_coverage, median_coverage, q1


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

viral_taxids = [
    10760,
    532076,
    1176767,
    1176765,
    1195074,
    1176434,
    227720,
    1176766,
    1837842,
    482822,
    1527506,
    2053563,
    1075775,
    2079317,
    1075774,
    10759,
    866889,
    1871708
]

sample_info = pd.read_excel('190904_Urine_Sample_Processing_Log.xlsx')
sample_info = sample_info.rename(columns={'Accession #': 'Accession'})
fqo = pd.read_csv('FastQataloguer_ARUP_Urine_2019-09-11.csv')
sample_info = sample_info.rename(columns={'Accession #': 'Accession'})
fqo['Accession'] = fqo['Accession'].str.upper()
fqo = fqo[~fqo['Run Directory/Run ID'].isna()]
quantifications = pd.read_csv('quantifications.csv')
merged = sample_info.merge(fqo[['Accession', 'Reads PostQual']], on='Accession')
merged = merged.merge(quantifications, on='Accession')
merged['taxid'] = merged['Detected Organism'].copy()
for organism in merged['Detected Organism'].unique():
    org_match_idx = merged.index[merged['Detected Organism'] == organism]
    for idx in org_match_idx:
        merged.at[merged.index[idx], 'taxid'] = org_taxids[organism]
summary_file_dir = '/Users/jmontgomery/Desktop/tmp_summary'
accession_ls, ctrl_ls, total_reads_ls, normalized_ctrl_ls, read_count_ls, \
    quantification_ls, coverage_mean_ls, coverage_median_ls, q1_ls = \
        [], [], [], [], [], [], [], [], []
for idx, row in merged.iterrows():
    accession = row['Accession']
    total_reads = row['Reads PostQual']
    quantification = 10**(row['log(Genomic Equivalents / ml)'])
    print(accession)
    # Get summary files
    summary_file_paths = glob.glob(os.path.join(summary_file_dir, '*' + accession + '*'))
    summary_file_paths = [path for path in summary_file_paths if 'dxsm.out.summary' in path]
    summary_file_paths = [path for path in summary_file_paths if not path.endswith('.done')]
    if len(summary_file_paths) < 3:  # Skip if summary files not present
        continue
    viral_path = [path for path in summary_file_paths if 'rna.viral' in path][0]
    viral_summary = read_summary_files(viral_path)
    bacterial_path = [path for path in summary_file_paths if 'rna.bacterial' in path][0]
    bacterial_summary = read_summary_files(bacterial_path)
    # Get ctrl counts
    ctrl_counts = []
    for org_info in viral_summary:
        if org_info['taxid'] in viral_taxids:
            ctrl_counts.append(org_info['read_count'])
    if len(ctrl_counts) < 1:
        print(f"No IC found for {accession}.")
        continue
    ctrl_counts_sum = np.sum(ctrl_counts)
    normalized_ctrl = 1e7 * (ctrl_counts_sum / total_reads)
    # Get quantifications
    read_count, mean_coverage, median_coverage, q1 = \
        get_read_counts(bacterial_summary, row['taxid'])
    accession_ls.append(accession)
    ctrl_ls.append(ctrl_counts_sum)
    normalized_ctrl_ls.append(normalized_ctrl)
    total_reads_ls.append(total_reads)
    read_count_ls.append(read_count)
    quantification_ls.append(quantification)
    coverage_median_ls.append(median_coverage)
    coverage_mean_ls.append(mean_coverage)
    q1_ls.append(q1)

df = pd.DataFrame(
    data={
        'Accession': accession_ls,
        'Ctrl Count': ctrl_ls,
        'Normalized IC': normalized_ctrl_ls,
        'Total Reads': total_reads_ls,
        'Organism Read Count': read_count_ls,
        'Quantification': quantification_ls,
        'Median Coverage': coverage_median_ls,
        'Mean Coverage': coverage_mean_ls,
        'Coverage Q1': q1_ls
    }
)

df['Detection'] = ''
df.loc[~df['Quantification'].isna(), 'Detection'] = 'TP'
fn_below_cutoff = ['IDBD-D100384', 'IDBD-D100413', 'IDBD-D100447', 'IDBD-D100450']
fn_no_detection = ['IDBD-D100414', 'IDBD-D100415', 'IDBD-D100416', 'IDBD-D100440', 'IDBD-D100441']
df.loc[df['Accession'].isin(fn_below_cutoff), 'Detection'] = 'FN - Below cutoff'
df.loc[df['Accession'].isin(fn_no_detection), 'Detection'] = 'FN - No detection'
df['log Quantification'] = np.log(df['Quantification'])
fig = px.box(df, y='Normalized IC', log_y=True, points='all',
             title='Normalized IC count, ARUP urine samples (n=96)',
             color='Detection', hover_data=['Accession', 'Detection', 'Normalized IC'])
fig.update_yaxes(range=[0, 5.1])
fig.show()
df.to_csv('Normalized_IC_Counts.csv')

