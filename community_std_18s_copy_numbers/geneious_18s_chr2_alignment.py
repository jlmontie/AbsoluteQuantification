import pandas as pd
import plotly.express as px

cov_df = pd.read_csv('data/CommunityStandard16sChr2Aln.csv')
cov_df['accession'] = cov_df['Name'].apply(lambda x: '-'.join(x.split('-')[3:5]))
copy_numbers = pd.DataFrame(cov_df.loc[cov_df['Reference'] == '18S', 'Mean Coverage'].values /
    cov_df.loc[cov_df['Reference'] == 'Chr2', 'Mean Coverage'].values, columns=['18S Copy Numbers'],
    index=cov_df['accession'].drop_duplicates())

fig = px.histogram(copy_numbers, x='18S Copy Numbers', nbins=500)
fig.show()