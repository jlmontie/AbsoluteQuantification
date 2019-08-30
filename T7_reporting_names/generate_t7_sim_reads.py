import subprocess as sp
from simulation_scripts.sim_reads import simulate_reads_tiled_by_depth as simulate
from Bio import SeqIO

seq_ls = []
for record in SeqIO.parse('T7_reporting_names/t7.fa', 'fasta'):
    seq_ls.append(str(record.seq))
sequence = ''.join(seq_ls)
len(sequence)

outfile_handle = open('T7_reporting_names/t7_simulated_reads.fa', 'w')
simulate('T7', sequence, 1000, 100, outfile_handle)