import pandas as pd
import plotly.express as px
import numpy as np
from scipy.stats import linregress


def get_trendline(df):
    lin_fit_x = df.loc[~(df['Plate count log(cfu)/ml'].isna()),
                          'Plate count log(cfu)/ml']
    lin_fit_y = df.loc[~(df['Plate count log(cfu)/ml'].isna()),
                          'log(Genomic Equivalents)/ml']
    slope, intercept, rval, pval, stderr = linregress(lin_fit_x, lin_fit_y)
    trendline_x=np.linspace(lin_fit_x.min(), lin_fit_x.max())
    trendline_y=slope * np.linspace(lin_fit_x.min(), lin_fit_x.max()) + intercept
    return trendline_x, trendline_y, rval, (slope, intercept)


ge_df = pd.read_csv('ge_cfu_correlation_plot_data.csv')
trendline_x, trendline_y, rval, coeff = get_trendline(ge_df)
# Validation report figure
fig1 = px.scatter(
    ge_df,
    x='Plate count log(cfu)/ml',
    y='log(Genomic Equivalents)/ml',
    # marginal_x='histogram', marginal_y='histogram',
    # color='Detected Organism',
    hover_data=['Accession',
                'Detected Organism',
                'Plate count log(cfu)/ml',
                'log(Genomic Equivalents)/ml'
    ]
)
fig1.add_scatter(x=np.linspace(4, 7), y=coeff[0] * np.linspace(4, 7) + coeff[1],
                     mode='lines', showlegend=False)
fig1.add_scatter(x=[4, 7], y=[7, 7],
                       mode='lines', line=dict(color='rgba(0.5, 0.5, 0.5, 0.8)',
                       width=4), showlegend=False)
fig1.add_scatter(x=[4, 7], y=[8, 8],
                       mode='lines', line=dict(color='rgba(0.5, 0.5, 0.5, 0.8)',
                       width=4), showlegend=False)
fig1.add_scatter(x=[5, 5], y=[4.5, 10.5],
                       mode='lines', line=dict(color='rgba(0.5, 0.5, 0.5, 0.8)',
                       width=4), showlegend=False)
fig1.add_scatter(x=[6, 6], y=[4.5, 10.5],
                       mode='lines', line=dict(color='rgba(0.5, 0.5, 0.5, 0.8)',
                       width=4), showlegend=False)
fig1.update_layout(
    yaxis=dict(
        range=[5, 10],
        tickvals = [5, 6, 7, 8, 9, 10],
        title='Explify v2.1 Platform Quantification [GE/mL]',
        ticktext=[
            '10<sup>5</sup>',
            '10<sup>6</sup>',
            '10<sup>7</sup>',
            '10<sup>8</sup>',
            '10<sup>9</sup>',
            '10<sup>10</sup>',
        ],
        tickfont=dict(size=24),
        titlefont=dict(size=24)
    ),
    xaxis=dict(
        title='Urine Culture Quantification [CFU/mL]',
        tickvals = [4, 5, 6, 7],
        ticktext=[
            '10<sup>4</sup>',
            '10<sup>5</sup>',
            '10<sup>6</sup>',
            '10<sup>7</sup>',
        ],
        tickfont=dict(size=24),
        titlefont=dict(size=24)
    ),
    height=850, width=600
)
fig1.write_html('Vertical.html')

# Horizontal
fig2 = px.scatter(
    ge_df,
    y='Plate count log(cfu)/ml',
    x='log(Genomic Equivalents)/ml',
    # marginal_x='histogram', marginal_y='histogram',
    # color='Detected Organism',
    hover_data=['Accession',
                'Detected Organism',
                'Plate count log(cfu)/ml',
                'log(Genomic Equivalents)/ml'
    ]
)
fig2.add_scatter(y=np.linspace(4, 7), x=coeff[0] * np.linspace(4, 7) + coeff[1],
                     mode='lines', showlegend=False)
fig2.add_scatter(y=[4, 7], x=[7, 7],
                       mode='lines', line=dict(color='rgba(0.5, 0.5, 0.5, 0.8)',
                       width=4), showlegend=False)
fig2.add_scatter(y=[4, 7], x=[8, 8],
                       mode='lines', line=dict(color='rgba(0.5, 0.5, 0.5, 0.8)',
                       width=4), showlegend=False)
fig2.add_scatter(y=[5, 5], x=[4.5, 10.5],
                       mode='lines', line=dict(color='rgba(0.5, 0.5, 0.5, 0.8)',
                       width=4), showlegend=False)
fig2.add_scatter(y=[6, 6], x=[4.5, 10.5],
                       mode='lines', line=dict(color='rgba(0.5, 0.5, 0.5, 0.8)',
                       width=4), showlegend=False)
fig2.update_layout(
    xaxis=dict(
        range=[5, 10],
        tickvals = [5, 6, 7, 8, 9, 10],
        title='Explify v2.1 Platform Quantification [GE/mL]',
        ticktext=[
            '10<sup>5</sup>',
            '10<sup>6</sup>',
            '10<sup>7</sup>',
            '10<sup>8</sup>',
            '10<sup>9</sup>',
            '10<sup>10</sup>',
        ],
        tickfont=dict(size=24),
        titlefont=dict(size=24)
    ),
    yaxis=dict(
        title='Urine Culture Quantification [CFU/mL]',
        range=[4, 7],
        tickvals = [4, 5, 6, 7],
        ticktext=[
            '10<sup>4</sup>',
            '10<sup>5</sup>',
            '10<sup>6</sup>',
            '10<sup>7</sup>',
        ],
        tickfont=dict(size=24),
        titlefont=dict(size=24)
    ),
    height=575, width=900
)
fig2.write_html('Horizontal.html')

# Square
fig3 = px.scatter(
    ge_df,
    x='Plate count log(cfu)/ml',
    y='log(Genomic Equivalents)/ml',
    # marginal_x='histogram', marginal_y='histogram',
    # color='Detected Organism',
    hover_data=['Accession',
                'Detected Organism',
                'Plate count log(cfu)/ml',
                'log(Genomic Equivalents)/ml'
    ]
)
fig3.add_scatter(x=np.linspace(4, 10), y=coeff[0] * np.linspace(4, 10) + coeff[1],
                     mode='lines', showlegend=False)
fig3.add_scatter(x=[4, 10], y=[7, 7],
                       mode='lines', line=dict(color='rgba(0.5, 0.5, 0.5, 0.8)',
                       width=4), showlegend=False)
fig3.add_scatter(x=[4, 10], y=[8, 8],
                       mode='lines', line=dict(color='rgba(0.5, 0.5, 0.5, 0.8)',
                       width=4), showlegend=False)
fig3.add_scatter(x=[5, 5], y=[4, 10],
                       mode='lines', line=dict(color='rgba(0.5, 0.5, 0.5, 0.8)',
                       width=4), showlegend=False)
fig3.add_scatter(x=[6, 6], y=[4, 10],
                       mode='lines', line=dict(color='rgba(0.5, 0.5, 0.5, 0.8)',
                       width=4), showlegend=False)
fig3.update_layout(
    yaxis=dict(
        range=[4, 10],
        tickvals = [4, 5, 6, 7, 8, 9, 10],
        title='Explify v2.1 Platform Quantification [GE/mL]',
        ticktext=[
            '10<sup>4</sup>',
            '10<sup>5</sup>',
            '10<sup>6</sup>',
            '10<sup>7</sup>',
            '10<sup>8</sup>',
            '10<sup>9</sup>',
            '10<sup>10</sup>',
        ],
        tickfont=dict(size=24),
        titlefont=dict(size=24)
    ),
    xaxis=dict(
        title='Urine Culture Quantification [CFU/mL]',
        range=[4, 10],
        tickvals = [4, 5, 6, 7, 8, 9, 10],
        ticktext=[
            '10<sup>4</sup>',
            '10<sup>5</sup>',
            '10<sup>6</sup>',
            '10<sup>7</sup>',
            '10<sup>8</sup>',
            '10<sup>9</sup>',
            '10<sup>10</sup>',
        ],
        tickfont=dict(size=24),
        titlefont=dict(size=24)
    ),
    height=850, width=850
)
fig3.write_html('Square.html')