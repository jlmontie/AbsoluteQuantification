import json
import glob
import os
import numpy as np
import pandas as pd

report_dir = 'synergy_validation_jan_2020'
classification_output_dir = '/data/analysis_group1/synergy_validation/quantification/validation_jan_2020/classification_pipeline_output'

path_ls = glob.glob(os.path.join(classification_output_dir, '*.rna.*.dxsm.out.summary'))
path_ls = [path for path in path_ls if not path.endswith('.done')]

taxid_ls = [562, 573, 5476, 5482]
sample_ls = []
org_ls = []
quant_ls = []
quant_log_ls = []
quant_bin_ls = []
coverage_ls = []
coverage_q1_ls = []
coverage_max_ls = []
coverage_min_ls = []
coverage_median_ls = []
coverage_q3_ls = []
for path in path_ls:
    with open(path) as file:
        for line in file:
            obj = json.loads(line)
            if obj['taxid'] in taxid_ls:
                sample_ls.append(os.path.basename(path).split('.')[0])
                org_ls.append(obj['name'])
                quant_ls.append(obj['absolute_quant'])
                if obj['absolute_quant'] is not None:
                    quant_log_ls.append(np.log10(obj['absolute_quant']))
                else:
                    quant_log_ls.append(None)
                quant_bin_ls.append(obj['absolute_quant_bin'])
                for gene in obj['gene_info']:
                    if gene['geneid'] == 0:
                        coverage_ls.append(gene['coverage'])
                        coverage_string = gene['coverage_string']
                        coverage = coverage_string.split(',')[:-1]
                        coverage = np.array([int(item) for item in coverage])
                        coverage_q1_ls.append(np.percentile(coverage, 25))
                        coverage_max_ls.append(np.max(coverage))
                        coverage_min_ls.append(np.min(coverage))
                        coverage_median_ls.append(np.percentile(coverage, 50))
                        coverage_q3_ls.append(np.percentile(coverage, 75))

df = pd.DataFrame(data={
    'sample': sample_ls,
    'organism': org_ls,
    'quantification': quant_ls,
    'quantification log': quant_log_ls,
    'quantification bin': quant_bin_ls,
    'coverage': coverage_ls,
    'depth min': coverage_min_ls,
    'depth q1': coverage_q1_ls,
    'depth median': coverage_median_ls,
    'depth q3': coverage_q3_ls,
    'depth max': coverage_max_ls
})
df.to_csv(os.path.join(report_dir, 'in_silico_quantification_results.csv'), index=False)
