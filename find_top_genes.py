import json
import numpy as np


def sort_genes(count_data, organism):
    gene_ls = []
    count_max_ls = []
    org_data = count_data[str(organism)]
    for gene in org_data["Gene Counts"].keys():
        gene_ls.append(gene)
        count_max_ls.append(np.array(org_data['Gene Counts'][gene]['Read Counts']).max())

    sort_idx = np.argsort(np.array(count_max_ls))
    sorted_genes = np.array(gene_ls)[sort_idx]
    return sorted_genes[::-1]