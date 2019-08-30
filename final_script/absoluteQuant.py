import json
import math

def absoluteQuant(ctrlReadCounts, summaryObjectList, rdnaCopyNumbers):
    slope = 1
    intercept = 8
    ctrlReadCountsMean = sum(ctrlReadCounts) / len(ctrlReadCounts)
    updatedSummaryObjectList = []
    for coverageInfo in summaryObjectList:
        taxid = str(coverageInfo['taxid'])
        if taxid not in rdnaCopyNumbers:  # skip if organism not in resource
            genomicEquivalents = 'NaN'
        else:
            for gene in coverageInfo['gene_info']:
                if gene['geneid'] == 0:
                    coverageStr = gene['coverage_string']
                    coverage = coverageStr.split(',')[:-1]
                    coverage = [int(item) for item in coverage]
                    coverage.sort()
                    coverageQ1 = coverage[:round((len(coverage) + 1) / 4)]  # Lower quartile of coverage
                    coverageQ1Mean = sum(coverageQ1) / len(coverageQ1)
                    rdnaAdjustment = rdnaCopyNumbers[taxid]['copies']
                    coverageNormalized = coverageQ1Mean / (rdnaAdjustment * ctrlReadCountsMean)
                    if coverageNormalized == 0:
                        genomicEquivalents = 0
                        continue
                    coverageLog = math.log10(coverageNormalized)  # Model is log-log
                    genomicEquivalentsLog = slope * coverageLog + intercept
                    genomicEquivalents = 10**genomicEquivalentsLog
        coverageInfo.update({'absolute_quant': genomicEquivalents})
        updatedSummaryObjectList.append(coverageInfo)
    return updatedSummaryObjectList
