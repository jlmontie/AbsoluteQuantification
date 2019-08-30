import json
from ncbi_taxonomy_utils import ncbi_taxonomy
import numpy as np
from collections import defaultdict
ncbi = ncbi_taxonomy()


def get_ctrl_counts(file_vir, ctrl_orgs):
    ctrl_count_ls = []
    for line in file_vir:
        cov_info_vir = json.loads(line)
        if cov_info_vir['name'].lower() in ctrl_orgs:
            ctrl_count_ls.append(cov_info_vir['read_count'])
    if len(ctrl_count_ls) > 0:
        ctrl_count = np.mean(ctrl_count_ls)
    else:
        ctrl_count = np.nan
    return ctrl_count


def get_organism_counts(file, target_taxids, ctrl_count=None):
    org_count_dict_ls = []
    for line in file:
        cov_info = json.loads(line)
        if cov_info['taxid'] in target_taxids:
            taxid = cov_info['taxid']
            organism = ncbi.get_name(taxid)
            gene_info = cov_info['gene_info']
            for gene in gene_info:
                if gene['geneid'] == 0:
                    read_count = gene['read_count']
            org_count_dict = {
                "taxid": taxid,
                "Organism": organism,
                "Read Counts": read_count,
                "Ctrl Counts": ctrl_count
            }
            if ctrl_count is not None:
                org_count_dict.update({"Ctrl Counts": ctrl_count})
            org_count_dict_ls.append(org_count_dict)
    return org_count_dict_ls


def combine_dictionaries(dict_ls):
    d = defaultdict(lambda: defaultdict(list))
    for dictionary in dict_ls:
        key_ls = list(dictionary.keys())
        key_ls.remove('taxid')
        for key in key_ls:
            d[dictionary['taxid']][key].append(dictionary[key])
    return d


def get_counts(file, target_taxids, file_vir=None, ctrl_orgs=None):
    if file_vir is not None:
        ctrl_count = get_ctrl_counts(file_vir, ctrl_orgs)
        org_count_dict_ls = get_organism_counts(file, target_taxids, ctrl_count)
    else:
        org_count_dict_ls = get_organism_counts(file, target_taxids)
    combined_count_dict = combine_dictionaries(org_count_dict_ls)
    return org_count_dict_ls
