import sys
sys.path.insert(0, '/uufs/chpc.utah.edu/common/home/u0002613/taxonomer/utils/explify_v2/simulation_scripts')
from sim_reads import simulate_reads_tiled_by_depth as simulate
import glob
import numpy as np

# Calculate required coverages
ic_count = 1e6
rdna_copies_dict = {
    '562_e_coli': 7.0,
    '573_k_pneumoniae': 8.0,
    '5476_c_albicans': 48.4,
    '5482_c_tropicalis': 48.4
}
# Model coefficients
m = 1.013875881594209
b = 12.446039813516776
quantification_ls = [10**(6.5), 10**(7.5), 10**(8.5)]

fa_ls = glob.glob("../ref_seqs/*.fa")
for fa in fa_ls:
    sequence_name = fa.split('/')[-1].split('.')[0]
    print(sequence_name)
    rdna_copies = rdna_copies_dict[sequence_name]
    with open(fa) as file_in:
        for quantification in  quantification_ls:
            depth = ic_count * rdna_copies * 10 ** ((np.log10(quantification) - b) / m)
            print(f"quantification: {quantification}")
            print(f"depth: {depth}")
            file_in.readline()
            sequence = file_in.readline()
            with open(f"{sequence_name}_quant_{quantification:0.0f}_depth_{depth:0.0f}", "w") as file_out:
                simulate(sequence_name, sequence, depth, 100, file_handle=file_out)