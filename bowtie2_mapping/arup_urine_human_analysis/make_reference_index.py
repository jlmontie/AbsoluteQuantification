import sys
sys.path.insert(0,'../')
import pipeline
import os
import subprocess as sp

for ref in os.listdir('reference_fastas'):
    fasta = os.path.join('reference_fastas', ref)
    pipeline.index_reference(fasta)
    sp.call(f"mv index_files index_{ref.split('.')[0]}",  shell=True)
