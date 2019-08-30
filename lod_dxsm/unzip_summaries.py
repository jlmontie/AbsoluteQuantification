import subprocess as sp
import os

parentdir = 'LOD_summary_files_dna'
summary_files = os.listdir(parentdir)
for file in summary_files:
    sp.Popen(['gunzip', os.path.join(parentdir, file)])