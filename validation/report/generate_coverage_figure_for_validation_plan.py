import json
import plotly.graph_objects as go
from plotly.colors import DEFAULT_PLOTLY_COLORS as colors
import numpy as np
import pandas as pd

# var_regions = 'V2::37::142::1,V3::333::397::1,V4::476::582::1,V5::722::779::1,V6::886::943::1,V7::1017::1073::1,V8::1143::1194::1'
var_regions = [
    [37, 142],
    [333, 397],
    [476, 582],
    [722, 779],
    [886, 943],
    [1017, 1073],
    [1143, 1194]
]
var_color = colors[2]
summary_path = '/Users/jmontgomery/OneDrive/Documents/IDbyDNA/Code/AbsoluteQuantification/synergy/zymo_2_titration_synergy/summary_with_quant/20190917SEQ-MCS2-1-d-35-0014-GGTAACTCGC-TCACCAACTT-5d851e6a-r.rna.bacterial.dxsm.out.summary'
with open(summary_path) as file:
    for idx, line in enumerate(file):
        obj = json.loads(line)
        if obj['taxid'] == 28901:
            cov_str = [int(cov) for cov in obj['gene_info'][0]['coverage_string'].split(',')[:-1]]
            cov_str[0], cov_str[-1] = 0, 0
            q1 = np.quantile(cov_str, 0.25)
            fig = go.Figure()
            for var in var_regions:
                fig.add_trace(
                    go.Scatter(
                        x=[var[0], var[0], var[1], var[1]],
                        y=[0, 16, 16, 0],
                        fill='toself',
                        mode='lines',
                        fillcolor='rgba(0.85, 0.85, 0.85, 0.5)',
                        line_color='rgba(0.85, 0.85, 0.85)',
                        showlegend=False
                    )
                )
            fig.add_scatter(y=cov_str, fill='tozerox', showlegend=False, fillcolor='rgb(79, 137, 185)', line_color=colors[0])  # 'rgba(31, 119, 180, 0.8)'
            print(colors[0])
            fig.add_scatter(x=[0, 1250], y=[q1, q1], showlegend=False, mode='lines', line_color=colors[1])
            fig.update_xaxes(title_text='Nucleotide position', title_font=dict(size=18), tickfont=dict(size=18), range=[0, 1250])
            fig.update_yaxes(title_text='Coverage depth', title_font=dict(size=18), tickfont=dict(size=18), range=[0, 15])
            fig.update_layout(
                title=obj['name'],
                annotations=[
                    dict(
                        x=790,
                        y=q1,
                        text='<b>Q1 coverage depth</b>',
                        arrowhead=2,
                        arrowsize=1,
                        arrowwidth=3,
                        ax=120,
                        ay=-135,
                        font=dict(size=24)
                    )
                ],
                height=600,
                width=1500,
                template='plotly_white'
            )
            fig.show()
            break
