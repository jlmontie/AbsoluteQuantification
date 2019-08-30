import pandas as pd
import plotly.express as px

# cfu_df = pd.read_csv('arup_organism_list/quantified_list.txt', sep='\t')
cfu_df = pd.read_csv('arup_organism_list/quantified_list_selected_single_positives.txt', sep='\t')
fig = px.histogram(cfu_df, x='concentration', labels={'concentration': 'concentration (cfu/ml)'},
                   histnorm='percent')
fig.show()

fig_org = px.histogram(cfu_df, x='organism', histnorm='percent').update_xaxes(categoryorder='sum descending')
fig_org.show()

detection_numbers = cfu_df.groupby('sample').size().to_frame()
detection_numbers = detection_numbers.rename(columns={0: 'number of detections'})
fig_num_detections = px.histogram(detection_numbers, x='number of detections', histnorm='percent')
fig_num_detections.show()