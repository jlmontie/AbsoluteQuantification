import scipy
import scikits.bootstrap as bootstrap
import pandas as pd

model_fit = pd.read_csv('synergy_second_std_curve/model_output_no_ecoli/plot_data.csv')
log_conc = model_fit['log_conc']
log_resid = model_fit['residuals']
resid_frac = 100 * (log_resid / log_conc)
# conc = 10**log_conc
# resid = 10**(log_conc - log_resid)
# resid.hist()
print(resid_frac)
resid_frac.hist()
CIs = bootstrap.ci(data=resid_frac, statfunction=scipy.median, alpha=0.05,
                   n_samples=100000)
print(CIs)
SEM = scipy.stats.sem(resid_frac)
print(SEM)