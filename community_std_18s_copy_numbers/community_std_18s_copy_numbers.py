from get_organism_counts import get_counts
import json

rna_dxsm_vir = open('lod_dxsm/community_std_titration/190725-1-1-idbd-d100291-d-07-ahmvfgafxy-cctctacatg-aggaggtatc-5d3b0b9c.rna.viral.dxsm.out.summary')
rna_dxsm_fungpar = open('lod_dxsm/community_std_titration/190725-1-1-idbd-d100291-d-07-ahmvfgafxy-cctctacatg-aggaggtatc-5d3b0b9c.rna.fungal_parasite.dxsm.out.summary')
dna_dxsm_vir = open('lod_dxsm/190725-1-1-IDBD-D100291-d-07-AHMVFGAFXY-CCTCTACATG-AGGAGGTATC-5d3b0b9c.dna.viral.dxsm.out.summary')
dna_dxsm_fungpar = open('lod_dxsm/190725-1-1-IDBD-D100291-d-07-AHMVFGAFXY-CCTCTACATG-AGGAGGTATC-5d3b0b9c.dna.fungal_parasite.fungal.dxsm.out.summary')

ctrl_orgs = ['enterobacteria phage t7', 'enterobacteria phage pr772']
target_taxids = [4932, 5207]
rna_counts = get_counts(rna_dxsm_fungpar, target_taxids, file_vir=rna_dxsm_vir, ctrl_orgs=ctrl_orgs)
dna_counts = get_counts(dna_dxsm_fungpar, target_taxids, file_vir=dna_dxsm_vir, ctrl_orgs=ctrl_orgs)