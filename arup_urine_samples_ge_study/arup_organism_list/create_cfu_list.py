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

sample_info_path = 'arup_urine_samples_ge_study/ge_distribution/190904_Urine_Sample_Processing_Log.xlsx'
samples = pd.read_excel(sample_info_path)
outpath = 'arup_urine_samples_ge_study/ge_distribution/sample_info.txt'
pattern = '>?\d{1,3}(,\d{3})\W\S+\W\S+\W\S+'
match = samples['RESULT LONG TEXT'].apply(search_string)
with open(outpath, 'w') as file:
    header = 'accession\tconcentration\torganism\tinfection_type\n'
    file.write(header)
    for (row_idx, row_data), (index, value) in zip(samples.iterrows(), match.iteritems()):
        for item in value:
            item_components = item.split(' ')
            concentration = item_components[0]
            concentration = concentration.replace('>', '')
            concentration = concentration.replace(',', '')
            accession = row_data['Accession #']
            organism = ' '.join(item_components[2:])
            if len(value) > 1:
                infection_type = 'multiple infection'
            else:
                infection_type = 'single infection'
            file.write(f"{accession}\t{concentration}\t{organism}\t{infection_type}\n")
