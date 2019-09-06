import json
import numpy as np
import plotly
import plotly.graph_objects as go
import gzip

# Used to build detection_disagreements.txt
# sample_summary = '/Users/jmontgomery/Desktop/tmp_summary_quant/190903-1-1-IDBD-D100414-d-04-AHMV3HAFXY-TCATAGATTG-CACCTTAATC-5d700753-r.rna.bacterial.dxsm.out.summary'
# idbd_detection_taxid = 287
# with open(sample_summary) as file:
#     for line in file:
#         obj = json.loads(line)
#         if obj['taxid'] == idbd_detection_taxid:
#             print(obj['absolute_quant'])
#             print(np.log10(obj['absolute_quant']))

# Get and plot coverage strings
sample_summary_ls = [
    '/uufs/chpc.utah.edu/common/home/idbydna-group3/results/idbd_dev/190903_NB551543_0129_AHMV3HAFXY/tax/190903-1-1-IDBD-D100414-d-04-AHMV3HAFXY-TCATAGATTG-CACCTTAATC-5d700753-r.rna.bacterial.dxsm.out.summary.gz',
    '/uufs/chpc.utah.edu/common/home/idbydna-group3/results/idbd_dev/190903_NB551543_0129_AHMV3HAFXY/tax/190903-1-1-IDBD-D100415-d-05-AHMV3HAFXY-GTATTCCACC-TTGTCTACAT-5d700753-r.rna.bacterial.dxsm.out.summary.gz',
    '/uufs/chpc.utah.edu/common/home/idbydna-group3/results/idbd_dev/190903_NB551543_0129_AHMV3HAFXY/tax/190903-1-1-IDBD-D100416-d-06-AHMV3HAFXY-CCTCCGTCCA-CACCGATGTG-5d700753-r.rna.bacterial.dxsm.out.summary.gz'
]
coverage_arr_ls = []
for summary_path in sample_summary_ls:
    with gzip.open(summary_path, 'rt') as file:
        for line in file:
            obj = json.loads(line)
            if obj['taxid'] == 287:
                gene_info = obj['gene_info']
                for gene in gene_info:
                    if gene['geneid'] == 0:
                        coverage_str = gene['coverage_string']
                        coverage_arr = [int(s) for s in coverage_str.split(',')[:-1]]
                        coverage_arr_ls.append(coverage_arr)
                    else:
                        print(summary_path)
                        print("16S not found")

accession_ls = ['IDBD-D100414', 'IDBD-D100415', 'IDBD-D100416']
data = []
for accession, coverage_arr in zip(accession_ls, coverage_arr_ls):
    trace = go.Scatter(
        y=coverage_arr,
        name=accession,
        mode='lines'
    )
    data.append(trace)
layout = go.Layout(
    xaxis=dict(
        title='Coverage'
    ),
    yaxis=dict(
        title='Position'
    )
)
fig = go.Figure(data, layout)
fig.write_html('paeruginosa_coverage.html')
