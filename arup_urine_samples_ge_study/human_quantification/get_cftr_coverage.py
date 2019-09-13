import pandas as pd
import sys
sys.path.insert(0, '../../bowtie_mapping')
import pipeline

pipeline.index_reference()

fqo = pd.read_csv('../ge_distribution/FastQataloguer_ARUP_Urine_2019-09-11.csv')