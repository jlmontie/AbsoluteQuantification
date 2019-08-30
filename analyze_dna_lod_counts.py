import json

with open('data/lod_counts_dna_raw.json') as file:
    dna_counts = json.load(file)

with open('lod_dxsm/LOD_summary_files_dna/190119-1-1-IDBD-V000235-d-02-AHLHVGBGX9-GGTTGCGAGG-TTGCTCTATT-5c4764b4.dna.bacterial.dxsm.out.summary') as file:
    kpneumo_sample_summary = json.loads(file.readline())

orgs_taxid = {
    'Hinfluenzae': 727,
    'Bpertussis': 520,
    'Kpneumoniae': 573,
    'Saureus': 1280,
    'Nfarcinica': 37329,
    'Nwallacei': 480035
}

kpneumo = dna_counts['573']
kpneumo_reads = kpneumo['Read Counts']

genes = []
for conc_reads in kpneumo_reads:
    conc_reads_genes = []
    for gene_dict in conc_reads:
        genes.append(list(gene_dict.keys())[0])
        conc_reads_genes.append(list(gene_dict.keys())[0])
    len(set(conc_reads_genes))
kpneumo_genes = set(genes)

for gene in kpneumo_genes:
    

# There appears to be a gene 100 in the dna_counts file for kpneumo that is not
# in the sample kpneumo summary file. Don't know why!
kpneumo_sample_summary_genes = []
for gene in kpneumo_sample_summary['gene_info']:
    kpneumo_sample_summary_genes.append(str(gene['geneid']))
len(kpneumo_sample_summary_genes)
len(set(kpneumo_sample_summary_genes))
set(kpneumo_sample_summary_genes).difference(kpneumo_genes)
kpneumo_genes.difference(set(kpneumo_sample_summary_genes))

