from sim_reads import simulate_reads_tiled_by_depth as simulate
from sim_reads import simulate_reads_tiled as simulate_tiled
import glob
import numpy as np
import os

# banks directory
output_dir = '/data/analysis_group1/synergy_validation/quantification/validation_jan_2020/simulated_sequences/'

# Calculate required coverages
ic_count = 1e6
rdna_copies_dict = {
    '562_e_coli': 7.0,
    '573_k_pneumoniae': 8.0,
    '5476_c_albicans': 48.4,
    '5482_c_tropicalis': 48.4
}
# Model coefficients
# Initial Synergy validation
# m = 1.013875881594209
# b = 12.446039813516776
# Second Synergy validation Jan 2020
m = 0.7797973340233221
b = 11.247277836330708

quantification_ls = [10**(6.5), 10**(7.5), 10**(8.5)]

fa_ls = glob.glob("../ref_seqs/*.fa")
for fa in fa_ls:
    sequence_name = fa.split('/')[-1].split('.')[0]
    print(sequence_name)
    rdna_copies = rdna_copies_dict[sequence_name]
    for quantification in  quantification_ls:
        with open(fa) as file_in:
            depth = ic_count * rdna_copies * 10 ** ((np.log10(quantification) - b) / m)
            print(f"quantification: {quantification}")
            print(f"depth: {depth}")
            file_in.readline()
            sequence = file_in.readline()
            filename_out = f"{sequence_name}_quant_{np.log10(quantification):0.1f}_depth_{depth:0.0f}.fa"
            with open(os.path.join(output_dir, filename_out), "w") as file_out:
                simulate(sequence_name, sequence, depth, 100, file_handle=file_out)

# T7 simulation
# ic = '../ref_seqs/T7 (090918RESP115-18250402467-d, NC_001604).fasta'
# T4 simulation
ic = '../ref_seqs/t4_viral_seqs_190405_40236.fasta'
with open(ic) as ic_file:
    ic_file.readline()
    ic_sequence = ic_file.readline()
    ic_sequence_name = os.path.splitext(os.path.basename(ic))[0]
    ic_blocks_ls = []
    while len(ic_blocks_ls) < 1e6:
        block = simulate_tiled(ic_sequence_name, ic_sequence, 1, 100)
        ic_blocks_ls.extend(block)
    ic_blocks_ls = ic_blocks_ls[:1000000]
    with open(os.path.join(output_dir, ic_sequence_name + '.fa'), 'w') as ic_out:
        for item in ic_blocks_ls:
            ic_out.write(item)

# # T7 simulation
# t7 = '../ref_seqs/T7 (090918RESP115-18250402467-d, NC_001604).fasta'
# with open(t7) as t7_file:
#     t7_file.readline()
#     t7_sequence = t7_file.readline()
#     t7_sequence_name = 'T7_090918RESP115-18250402467-d_NC_001604'
#     t7_blocks_ls = []
#     while len(t7_blocks_ls) < 1e6:
#         block = simulate_tiled(t7_sequence_name, t7_sequence, 1, 100)
#         t7_blocks_ls.extend(block)
#     t7_blocks_ls = t7_blocks_ls[:1000000]
#     with open(os.path.join(output_dir, 'T7_090918RESP115-18250402467-d_NC_001604.fa'), 'w') as t7_out:
#         for item in t7_blocks_ls:
#             t7_out.write(item)