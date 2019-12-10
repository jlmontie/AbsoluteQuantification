import os
import json
from collections import defaultdict
import pandas as pd

orgs = {
    'IDBD-D100412': {
        'orgs': ['Escherichia coli'],
        'taxids': [562]
    },
    'IDBD-D100377': {
        'orgs': ['Escherichia coli'],
        'taxids': [562]
    },
    'IDBD-D100456': {
        'orgs': [
            'Escherichia coli',
            'Klebsiella pneumoniae',
            'Staphylococcus aureus'
        ],
        'taxids': [
            562,
            573,
            1280
        ]
    },
    'IDBD-D100482': {
        'orgs': ['Escherichia coli'],
        'taxids': [562]
    },
    'IDBD-D100395': {
        'orgs': [
            'Escherichia coli',
            'Proteus mirabilis',
            'Staphylococcus aureus'
        ],
        'taxids': [
            562,
            584,
            1280
        ]
    },
    'IDBD-D100382': {
        'orgs': [
            'Escherichia coli',
            'Staphylococcus aureus'
        ],
        'taxids': [
            562,
            1280
        ]
    },
    'IDBD-D100477': {
        'orgs': [
            'Escherichia coli',
            'Staphylococcus aureus'
        ],
        'taxids': [
            562,
            1280
        ]
    },
    'IDBD-D100438': {
        'orgs': [
            'Escherichia coli',
            'Klebsiella pneumoniae'
        ],
        'taxids': [
            562,
            573
        ]
    },
    'IDBD-D100461': {
        'orgs': ['Escherichia coli'],
        'taxids': [562]
    },
    'IDBD-D100400': {
        'orgs': ['Escherichia coli'],
        'taxids': [562]
    },
    '20190816-T7V1-C25-3-d-34': {
        'orgs': [
            'Pseudomonas aeruginosa',
            'Gardnerella vaginalis',
        ],
        'taxids': [
            287,
            2702
        ]
    },
    '20190820SEQ-C10-d-32': {
        'orgs': [
            'Staphylococcus aureus',
            'Gardnerella vaginalis',
            'Staphylococcus argenteus'
        ],
        'taxids': [
            1280,
            2702,
            985002
        ]
    },
    '20190813S-C27-1S1-d-12': {
        'orgs': [
            'Klebsiella pneumoniae',
            'Escherichia coli',
            'Gardnerella vaginalis'
        ],
        'taxids': [
            573,
            562,
            2702
        ]
    }
}

summary_parentdir = '/home/jmontgomery/AbsoluteQuantification/synergy/reproducibility/summary_with_quant20190820SEQ-C10-d-32'
output = '/home/jmontgomery/tmp/arup_urine_quants.csv'
files = os.listdir(summary_parentdir)
files = [file for file in files if any(key in file for key in list(orgs.keys()))]
print(files)
quants = defaultdict(lambda: defaultdict(list))
accession_ls = []
coverage_ls = []
organism_ls = []
quant_ls = []
for file in files:
    accession = [key for key in list(orgs.keys()) if key in file][0]
    taxids = orgs[accession]['taxids']
    with open(os.path.join(summary_parentdir, file)) as input:
        for line in input:
            org_info = json.loads(line)
            if org_info['taxid'] in taxids:
                quant = org_info['absolute_quant']
                for gene in org_info['gene_info']:
                    if gene['geneid'] == 0:
                        coverage = gene['coverage']
                quants[accession]['absolute_quant'].append(quant)
                quants[accession]['coverage'].append(coverage)
                quants[accession]['organism'].append(org_info['name'])
                accession_ls.append(accession)
                coverage_ls.append(coverage)
                organism_ls.append(org_info['name'])
                quant_ls.append(quant)


print(quants)
quants_df = pd.DataFrame(data={
    'Accession': accession_ls,
    'Organism': organism_ls,
    'Quantification': quant_ls,
    'Coverage': coverage_ls
})
quants_df.to_csv(output)
