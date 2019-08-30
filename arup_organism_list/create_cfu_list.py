import pandas as pd
import re

def search_string(str):
    pattern = '(>)?\d{1,3}(,\d{3})\W\S+\W\S+\W\S+'
    match_ls = []
    start_pos = 0
    for meaningless_idx in range(5):
        match_obj = re.search(pattern, str[start_pos:])
        if match_obj is None:
            break
        match = match_obj.group()
        match_ls.append(match)
        match_end = match_obj.span()[1]
        start_pos = match_end
    match_ls = list(set(match_ls))
    return match_ls

# samples = pd.read_excel('arup_organism_list/2019_08_26_UTI_Samples.xlsx',
#                         sheet_name='Samples')
samples = pd.read_csv('arup_organism_list/selected_samples_single_positives.csv')
pattern = '>?\d{1,3}(,\d{3})\W\S+\W\S+\W\S+'
match = samples['RESULT LONG TEXT'].apply(search_string)
# with open('arup_organism_list/quantified_list.txt', 'w') as file:
with open('arup_organism_list/quantified_list_selected_single_positives.txt', 'w') as file:
    header = 'sample\tconcentration\torganism\tinfection_type\n'
    file.write(header)
    for index, value in match.iteritems():
        for item in value:
            item_components = item.split(' ')
            concentration = item_components[0]
            concentration = concentration.replace('>', '')
            concentration = concentration.replace(',', '')
            organism = ' '.join(item_components[2:])
            if len(value) > 1:
                infection_type = 'multiple infection'
            else:
                infection_type = 'single infection'
            file.write(f"{index}\t{concentration}\t{organism}\t{infection_type}\n")
