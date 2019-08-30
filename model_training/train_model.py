import pandas as pd

fqo = pd.read_csv('model_training/FastQataloguer_CommunityStdTitration_2019-07-26.csv')
classification_dir = fqo['Diagnostic Output Dir']