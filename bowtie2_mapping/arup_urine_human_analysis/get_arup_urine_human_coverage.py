import sys
import pandas as pd
sys.path.insert(0,'../')
import pipeline

index_path = 'index_files'
fqo = pd.read_csv('../ge_distribution/FastQataloguer_ARUP_Urine_2019-09-11.csv')
