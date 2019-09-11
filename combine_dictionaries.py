from collections import defaultdict

def combine_dictionaries(dict_ls):
    d = defaultdict(lambda: defaultdict(list))
    for dictionary in dict_ls:
        key_ls = list(dictionary.keys())
        key_ls.remove('taxid')
        for key in key_ls:
            d[dictionary['taxid']][key].append(dictionary[key])
    return d