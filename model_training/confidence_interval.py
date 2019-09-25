import scipy
import scikits.bootstrap as bootstrap
import pandas as pd

model_fit = pd.read_csv('synergy_second_std_curve/model_output_no_ecoli/plot_data.csv')
log_conc = model_fit['log_conc']
log_resid = model_fit['residuals']
conc = 10**log_conc
resid = 10**(log_conc - log_resid)
resid.hist()
