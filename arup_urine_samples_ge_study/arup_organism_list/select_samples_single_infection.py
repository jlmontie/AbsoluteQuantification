import pandas as pd

sample_list = pd.read_csv('arup_organism_list/quantified_list.txt', sep='\t')
single_below_100k = sample_list[(sample_list['concentration'] < 1e5) &
                                (sample_list['infection_type'] == 'single infection')]
ecoli_single_100k = sample_list[(sample_list['concentration'] == 100000) &
                                (sample_list['organism'] == 'Escherichia coli') &
                                (sample_list['infection_type'] == 'single infection')]
sample_indexes = single_below_100k['sample'].tolist() + ecoli_single_100k['sample'].tolist()

sample_df = pd.read_excel('arup_organism_list/2019_08_26_UTI_Samples.xlsx',
                          sheet_name='Samples')
selected_samples = sample_df.iloc[sample_df.index[sample_indexes]]
selected_samples.to_csv('arup_organism_list/selected_samples_single_positives.csv')