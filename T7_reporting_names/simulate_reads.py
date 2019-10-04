import sys
sys.path.insert(0, '/uufs/chpc.utah.edu/common/home/u0002613/taxonomer/utils/explify_v2/simulation_scripts')
from sim_reads import simulate_reads_tiled as simulate_tiled
import glob
import numpy as np

t7 = '../validation/ref_seqs/T7 (090918RESP115-18250402467-d, NC_001604).fasta'
with open(t7) as t7_file:
    t7_file.readline()
    t7_sequence = t7_file.readline()
    t7_sequence_name = 'T7_090918RESP115-18250402467-d_NC_001604'
    t7_blocks_ls = []
    while len(t7_blocks_ls) < 1e6:
        block = simulate_tiled(t7_sequence_name, t7_sequence, 1, 100)
        t7_blocks_ls.extend(block)
    t7_blocks_ls = t7_blocks_ls[:1000000]
    with open('../ref_seqs/T7_090918RESP115-18250402467-d_NC_001604.fa', 'w') as t7_out:
        for item in t7_blocks_ls:
            t7_out.write(item)
