import pandas as pd
import numpy as np

cfu_df = pd.read_csv('arup_organism_list/quantified_list.txt', sep='\t')
sample_df = pd.read_excel('arup_organism_list/2019_08_26_UTI_Samples.xlsx',
                          sheet_name='Samples')
bin_vals = [1e3, 1e4, 2e4, 3e4, 4e4, 5e4, 6e4, 7e4, 8e4, 9e4, 10e5]
bin_num, bin_vals = np.histogram(cfu_df['concentration'], bins=bin_vals)
selection_num = 100 * bin_num / len(cfu_df)
selection_num = np.ceil(selection_num).astype(int)
assigned_bins = np.digitize(cfu_df['concentration'], bin_vals)
selection_df_ls = []
total_samples = 0
for idx, selection in enumerate(selection_num):
    bin_selection = idx + 1
    all_indices_at_bin = np.argwhere(assigned_bins == bin_selection).flatten()
    samples_at_indices = list(set(cfu_df.loc[cfu_df.index[all_indices_at_bin], 'sample'].tolist()))
    total_samples += selection
    if len(samples_at_indices) > 0:
        np.random.seed(999)
        if total_samples > 100:
            print(selection)
            selection = selection - (total_samples - 100) # ensure only 100 samples final
            print(selection)
        random_selection_of_indices = np.random.choice(samples_at_indices, selection, replace=False)
        selection_df_ls.append(sample_df.iloc[sample_df.index[random_selection_of_indices]])
pd.concat(selection_df_ls).to_csv('arup_organism_list/selected_samples.csv')
