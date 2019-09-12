import plotly.express as px
import pandas as pd
import numpy as np

df = pd.read_csv('~/Downloads/test.csv')
df = df.rename(columns={'Normalized Ctrl': 'Normalized IC'})
df['log Quantification'] = np.log(df['Quantification'])
# df['Detection'] = ''
# df.loc[df['Quantification'].isna(), 'Detection'] = 'False Negative'
# df.loc[~df['Quantification'].isna(), 'Detection'] = 'True Positive'
fig = px.box(df, y='Normalized IC', log_y=True, points='all',
             title='Normalized IC count, ARUP urine samples (n=96)',
             color='Detection', hover_data=['Accession', 'Detection', 'Normalized IC'])
fig.update_traces(marker_size=df['log Quantification'].tolist())
fig.update_yaxes(range=[0, 5.1])
fig.show()
