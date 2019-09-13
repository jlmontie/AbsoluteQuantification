import plotly.graph_objects as go
import plotly.express as px
import os
import json
import pandas as pd
import numpy as np

quant_files = os.listdir('simplified_quant')
study_info = pd.read_excel('20190812_Summary_of_Synergy_Reproducibility_Samples_Runs.xlsx')
batch_key = dict(zip(study_info['Synergy ID'].tolist(), study_info['Batch ID'].tolist()))
sample = []
batch = []
taxid = []
name = []
quant = []
for file in quant_files:
    path = os.path.join('simplified_quant', file)
    sample_id = '-'.join(file.split('-')[1:3])
    display_id = file.split('-')[1]
    sample.append(display_id)
    batch.append(batch_key[sample_id])
    if os.stat(path).st_size == 0:
        taxid.append(np.nan)
        name.append(np.nan)
        quant.append(np.nan)
        continue
    with open(path) as infile:
        for line in (infile):
            quant_dict = json.loads(line)
            taxid.append(quant_dict['taxid'])
            name.append(quant_dict['name'])
            quant.append(quant_dict['quant'])
df = pd.DataFrame(data={'Sample ID':sample, 'Batch': batch, 'taxid': taxid,
                  'Organism': name, 'quant': quant})
df['log(GE)/ml'] = np.log10(df['quant'])
df = df[~(df['Sample ID'].isin(['C14', 'C16']))]

df_inter = df[df['Sample ID'].isin(['C17', 'C18', 'C21', 'C22', 'C24', 'C25', 'C26', 'C27'])]
print(len(df_inter['quant'].isna()))
# fig_inter = px.scatter(df_inter, x='Sample ID', y='log(GE)/ml', color='Batch',
#                 category_orders={'Sample ID': df_inter['Sample ID'].sort_values()},
#                 hover_data=['Batch', 'Organism', 'log(GE)/ml'])
# fig_inter.show()

df_inter_norm = df_inter.copy()
inter_mean = df_inter_norm.groupby('Sample ID').transform('mean')
df_inter_norm['Deviation from sample mean, log(GE)/ml'] = \
    df_inter_norm['log(GE)/ml'] - inter_mean['log(GE)/ml']
fig_inter_norm = px.box(
    df_inter_norm, y='Deviation from sample mean, log(GE)/ml',
    x='Organism',
    hover_data=['Batch', 'Organism', 'log(GE)/ml'],
    points='all',
    category_orders={
        'Organism': [
            'Candida tropicalis',
            'Pseudomonas aeruginosa',
            'Cryptococcus laurentii (Papiliotrema laurentii)',
            'Escherichia coli',
            'Candida parapsilosis',
            'Enterococcus faecalis',
            'Candida albicans',
            'Klebsiella pneumoniae'
        ]
    }
)
fig_inter_norm.update_layout(title='Inter Run Variation')
fig_inter_norm.show()
fig_inter_norm.write_html('InterRunVariation.html')

fig_inter = px.scatter(df_inter[df['Organism'].isin(['Candida tropicalis',
                                                    'Pseudomonas aeruginosa'])],
                       x='Organism', y='log(GE)/ml', color='Batch',
                       hover_data=['Batch', 'Organism', 'log(GE)/ml'])
fig_inter.update_traces(marker=dict(size=14))
fig_inter.update_xaxes(range=[-1, 2])
fig_inter.update_layout(height=700, width=700)
fig_inter.show()
fig_inter.write_html('LargeInterRunVariation.html')

df_intra = df[df['Batch']=='20190813B']
intra_mean = df_intra.groupby('Sample ID').transform('mean')
df_intra['Deviation from sample mean, log(GE)/ml'] = \
    df_intra['log(GE)/ml'] - intra_mean['log(GE)/ml']
fig_intra = px.box(df_intra, y='Deviation from sample mean, log(GE)/ml',
                   hover_data=['Batch', 'Organism', 'log(GE)/ml'],
                   points='all')
fig_intra.update_layout(title='Intra Run Variation', height=700, width=700)
fig_intra.show()
fig_intra.write_html('IntraRunVariation.html')