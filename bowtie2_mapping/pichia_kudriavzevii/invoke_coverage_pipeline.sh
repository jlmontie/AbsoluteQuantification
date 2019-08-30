#!/usr/bin/env bash
pipeline_path=/Users/jmontgomery/OneDrive/Documents/IDbyDNA/Code/AbsoluteQuantification/bowtie2_mapping/coverage_calc_pipeline.py
ref_gen=/Users/jmontgomery/OneDrive/Documents/IDbyDNA/Code/AbsoluteQuantification/data/ref_genome/4909_GCF_003054445.1_ASM305444v1_genomic.fna.gz
ref_18s=/Users/jmontgomery/OneDrive/Documents/IDbyDNA/Code/AbsoluteQuantification/data/18s_fastas/4909.fa
fastq=/Users/jmontgomery/Downloads/20190816-T7V1-C15-1S1-d-02-0008-GGTTGCGAGG-TTGCTCTATT-5d59b4d3_S41_R1_001.fastq
python3 $pipeline_path -r $ref_18s -g $ref_gen -s $fastq