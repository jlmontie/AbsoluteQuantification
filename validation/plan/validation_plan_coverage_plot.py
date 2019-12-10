import json
import plotly.graph_objects as go
import numpy as np

summary_path = '../../synergy/zymo_2_titration_idbd/summary_with_quant/190917-1-1-IDBD-D100520-d-01-AHTCM2AFXY-GAACTGAGCG-CGCTCCACGA-5d827495-r.rna.bacterial.dxsm.out.summary'

with open(summary_path) as file:
    for line in file:
        obj = json.loads(line)
        if obj['taxid'] == 287: # 1639:
            for gene in obj['gene_info']:
                if gene['geneid'] == 0:
                    cov_str = gene['coverage_string']
coverage = cov_str.split(',')[:-1]
coverage = np.array([int(item) for item in coverage])
coverage[0] = 0
coverage[-1] = 0
q1 = np.percentile(coverage, 25)
fig = go.Figure()
fig.add_trace(go.Scatter(y=coverage, mode='lines', fill='tonext'))
fig.add_trace(go.Scatter(y=[q1, q1], x=[0, len(coverage)], mode='lines'))
fig.update_yaxes(range=[0, 75], title='Coverage Depth')
fig.update_xaxes(title='Nucleotide Position')
fig.show()
fig2 = go.Figure()
fig2.add_box(y=coverage)
fig2.show()
