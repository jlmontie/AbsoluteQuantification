import os
from skimage import io
from skimage.viewer import ImageViewer
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
from scipy.stats import linregress
from plotly.colors import DEFAULT_PLOTLY_COLORS as colors


def test(x):
    num_streaks = x['Streaks']
    new_cfu_ml = streak_len_adjustment.loc[
        streak_len_adjustment.index[:num_streaks + 1], 'colonies'].sum() * 1000
    return new_cfu_ml


manual_counts = pd.read_excel('wasp_image_counting.xlsx')
streak_len_adjustment = pd.read_excel('StreakLengthAdjustment.xlsx')
manual_counts['Adjusted plate count cfu/ml'] = manual_counts['cfu/ml'].copy()
manual_counts.loc[~manual_counts['Streaks'].isna(), 'Adjusted plate count cfu/ml'] = \
    manual_counts[~manual_counts['Streaks'].isna()].apply(test, axis=1)
manual_counts['Adjusted plate count log(cfu)/ml'] = \
    manual_counts['Adjusted plate count cfu/ml'].copy()
manual_counts.loc[~(manual_counts['Adjusted plate count cfu/ml'] == 0),
                  'Adjusted plate count log(cfu)/ml'] = \
    np.log10(manual_counts.loc[~(manual_counts['Adjusted plate count cfu/ml'] == 0),
                  'Adjusted plate count log(cfu)/ml'])
manual_counts.loc[(manual_counts['Adjusted plate count cfu/ml'] == 0),
                  'Adjusted plate count log(cfu)/ml'] = np.nan

ge_counts = pd.read_csv('../arup_urine_samples_ge_study/ge_distribution/quantifications.csv')
merged = manual_counts.merge(ge_counts, on='IDBD #')
merged_no_na = merged[~(merged['Adjusted plate count log(cfu)/ml'].isna()) &
                      ~(merged['log(Genomic Equivalents)/ml'].isna())]
slope_adj, intercept_adj, rval_adj, _, _ = \
    linregress(merged_no_na['Adjusted plate count log(cfu)/ml'],
               merged_no_na['log(Genomic Equivalents)/ml'])
fit_x_adj = np.linspace(merged_no_na['Adjusted plate count log(cfu)/ml'].min(),
                    merged_no_na['Adjusted plate count log(cfu)/ml'].max())
fit_y_adj = slope_adj * fit_x_adj + intercept_adj
slope, intercept, rval, _, _ = \
    linregress(merged_no_na['Plate count log(cfu)/ml'],
               merged_no_na['log(Genomic Equivalents)/ml'])
fit_x = np.linspace(merged_no_na['Plate count log(cfu)/ml'].min(),
                        merged_no_na['Plate count log(cfu)/ml'].max())
fit_y = slope * fit_x + intercept
fig_len_adj = go.Figure()
fig_len_adj.add_scatter(x=merged_no_na['Adjusted plate count log(cfu)/ml'],
                        y=merged_no_na['log(Genomic Equivalents)/ml'],
                        name='Streak length adjustment', mode='markers',
                        marker=dict(color=colors[0], size=10))
fig_len_adj.add_scatter(x=fit_x_adj, y=fit_y_adj, mode='lines',
                        marker=dict(color=colors[0]), showlegend=False)
fig_len_adj.add_scatter(x=merged_no_na['Plate count log(cfu)/ml'],
                        y=merged_no_na['log(Genomic Equivalents)/ml'],
                        name='No adjustment', mode='markers',
                        marker=dict(color=colors[1]))
fig_len_adj.add_scatter(x=fit_x, y=fit_y, mode='lines',
                        marker=dict(color=colors[1]), showlegend=False)
fig_len_adj.update_layout(
    annotations=[
        dict(x=4.75, y=7.5, text=f"Pearson r after adjustment:<br>{rval_adj:0.2f}",
             showarrow=False, font=dict(size=16)),
        dict(x=4.75, y=8, text=f"Pearson r before adjustment:<br>{rval:0.2f}",
             showarrow=False, font=dict(size=16))]
)
fig_len_adj.show()
test = px.get_trendline_results(fig_len_adj)
image_dir = 'cropped_images'
image_files = os.listdir(image_dir)
image_files = [file for file in image_files if file.endswith('png') or file.endswith('jpg')]
print(image_files[0])
img = io.imread(os.path.join(image_dir, image_files[0]))
print(np.sum(np.where(img < 60)))
idbd_num_ls = []
mean_green_blue_ls = []
with open('green_blue_pixels.txt', 'w') as outfile:
    outfile.write(f"sample\tg_b_pixel_avg\n")
    for file in image_files:
        idbd_num_ls.append(file.split('.')[0])
        img = io.imread(os.path.join(image_dir, file))
        red = img[:, :, 0].ravel()
        red = red[red < 250]
        # not_red = np.argwhere(red[red < 5e3])
        green = img[:, :, 1].ravel()
        green[(green > 60) & (green < 10)] = 0
        green_img = ImageViewer(green.reshape(img.shape[:2]))
        green_img.show()
        blue = img[:, :, 2].ravel()
        blue[(blue > 60) & (blue < 10)] = 0
        blue_img = ImageViewer(blue.reshape(img.shape[:2]))
        blue_img.show()
        # break
        mean_green_blue = np.mean([len(green), len(blue)])
        mean_green_blue_ls.append(mean_green_blue)
        # df = pd.DataFrame(data={'red': red, 'green': green, 'blue': blue})
        # fig = px.histogram(df, x='red')
        # fig.add_histogram(x=df['green'], name='green')
        # fig.add_histogram(x=df['blue'], name='blue')
        # fig.show()
        # print(img.shape)
        outfile.write(f"{file}\t{mean_green_blue}\n")
df_pixel = pd.DataFrame(data={'IDBD #': idbd_num_ls, 'mean pixels': mean_green_blue_ls})
df = manual_counts.merge(df_pixel)
df_streaks = df[~df['Streaks'].isna()]
fig = px.scatter(df_streaks, x='Streaks', y='mean pixels', title='Green blue channel pixels and manual counting')
fig.show()