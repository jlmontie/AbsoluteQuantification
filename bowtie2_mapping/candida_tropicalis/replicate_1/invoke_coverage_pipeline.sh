#!/usr/bin/env bash
pipeline_path=/Users/jmontgomery/OneDrive/Documents/IDbyDNA/Code/AbsoluteQuantification/bowtie2_mapping/coverage_calc_pipeline.py
ref_gen=/Users/jmontgomery/OneDrive/Documents/IDbyDNA/Code/AbsoluteQuantification/data/ref_genome/5482_GCF_000006335.3_ASM633v3_genomic.fna.gz
ref_18s=/Users/jmontgomery/OneDrive/Documents/IDbyDNA/Code/AbsoluteQuantification/data/18s_fastas/5482.fa
fastq=/Users/jmontgomery/Downloads/20190816-T7V1-C18-1S1-d-05-0008-GCCGCACTCT-CGAGGTCGGA-5d59b4d3_S44_R1_001.fastq.gz
python3 $pipeline_path -r $ref_18s -g $ref_gen -s $fastq