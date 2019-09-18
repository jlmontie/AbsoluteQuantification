import json
import numpy as np
from ncbi_taxonomy_utils import ncbi_taxonomy
ncbi = ncbi_taxonomy()


def get_species_genus_taxid(taxid):
    taxid = int(taxid)
    path = ncbi.get_path(taxid)
    ranks = [ncbi.get_rank(path_taxid) for path_taxid in path]
    if 'species' in ranks:
        species_idx = ranks.index('species')
        species_taxid = str(path[species_idx])
    else:
        species_taxid = None
    if 'genus' in ranks:
        genus_idx = ranks.index('genus')
        genus_taxid = str(path[genus_idx])
    else:
        genus_taxid = None
    return species_taxid, genus_taxid


def absoluteQuant(ctrlReadCounts, summaryObjectList, rdnaCopyNumbers, quant_mode='coverage'):
    if quant_mode == 'coverage':
        # slope = 0.9752025158801086
        # intercept = 9.94632629381537
        # 5 log urine extraction, no P aeruginosa
        slope = 1.2169746581769076
        intercept = 10.99130819075618
    elif quant_mode == 'read_count':
        slope = 0.9777720637374567
        intercept = 8.752073402926564

    ctrlReadCountsMean = sum(ctrlReadCounts) #/ len(ctrlReadCounts)
    copyNumberMean = np.mean([val['copies'] for val in rdnaCopyNumbers.values()])
    updatedSummaryObjectList = []
    for coverageInfo in summaryObjectList:
        taxid = str(coverageInfo['taxid'])
        if taxid in rdnaCopyNumbers:
            rdnaAdjustment = rdnaCopyNumbers[taxid]['copies']
        else:
            species_taxid, genus_taxid = get_species_genus_taxid(taxid)
            if species_taxid in rdnaCopyNumbers:
                rdnaAdjustment = rdnaCopyNumbers[species_taxid]['copies']
            elif genus_taxid in rdnaCopyNumbers:
                rdnaAdjustment = rdnaCopyNumbers[genus_taxid]['copies']
            else:
                rdnaAdjustment = copyNumberMean
        for gene in coverageInfo['gene_info']:
            if gene['geneid'] == 0:
                # if gene['coverage'] < 0.97:
                    # genomicEquivalents = np.nan
                    # continue
                if quant_mode == 'coverage':
                    coverageStr = gene['coverage_string']
                    coverage = coverageStr.split(',')[:-1]
                    coverage = np.array([int(item) for item in coverage])
                    # coverage = coverage[coverage > 0]  # this overestimates quant w/ lower coverage!
                    q1 = np.quantile(coverage, 0.25)
                    coverageQ1 = coverage[coverage <= q1]
                    coverageQ1Mean = np.mean(coverageQ1)
                    coverageNormalized = coverageQ1Mean / (rdnaAdjustment * ctrlReadCountsMean)
                    if coverageNormalized == 0:
                        genomicEquivalents = 0
                        continue
                    coverageLog = np.log10(coverageNormalized)  # Model is log-log
                    genomicEquivalentsLog = slope * coverageLog + intercept
                    genomicEquivalents = 10**genomicEquivalentsLog
                    break
                elif quant_mode == 'read_count':
                    read_count = gene['read_count']
                    if read_count == 0:
                        genomicEquivalents = 0
                        continue
                    read_count_normalized = read_count / (rdnaAdjustment * ctrlReadCountsMean)
                    read_count_log = np.log10(read_count_normalized)
                    genomicEquivalentsLog = slope * read_count_log + intercept
                    genomicEquivalents = 10**genomicEquivalentsLog
        coverageInfo.update({'absolute_quant': genomicEquivalents})
        updatedSummaryObjectList.append(coverageInfo)
    return updatedSummaryObjectList
