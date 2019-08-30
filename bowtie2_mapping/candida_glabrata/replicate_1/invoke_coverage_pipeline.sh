#!/usr/bin/env bash
pipeline_path=/Users/jmontgomery/OneDrive/Documents/IDbyDNA/Code/AbsoluteQuantification/bowtie2_mapping/coverage_calc_pipeline.py
ref_gen=/Users/jmontgomery/OneDrive/Documents/IDbyDNA/Code/AbsoluteQuantification/data/ref_genome/5478_GCF_000002545.3_ASM254v2_genomic.fna.gz
ref_18s=/Users/jmontgomery/OneDrive/Documents/IDbyDNA/Code/AbsoluteQuantification/data/18s_fastas/5478.fa
fastq=/Users/jmontgomery/Downloads/20190816-T7V1-C16-1S1-d-03-0008-TAAGCATCCA-AATGGATTGA-5d59b4d3_S42_R1_001.fastq.gz
python3 $pipeline_path -r $ref_18s -g $ref_gen -s $fastq