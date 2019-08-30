pipeline_path=/Users/jmontgomery/OneDrive/Documents/IDbyDNA/Code/AbsoluteQuantification/bowtie2_mapping/coverage_calc_pipeline.py
ref_gen=/Users/jmontgomery/OneDrive/Documents/IDbyDNA/Code/AbsoluteQuantification/data/ref_genome/GCF_000149245.1_CNA3_genomic.fna
ref_18s=/Users/jmontgomery/OneDrive/Documents/IDbyDNA/Code/AbsoluteQuantification/data/18s_fastas/5207.fa
fastq=/Users/jmontgomery/OneDrive/Documents/IDbyDNA/Code/AbsoluteQuantification/data/fastq/190812-1-1-IDBD-D100331-d-01-AHMV5GAFXY-CGGTTACGGC-CTATAGTCTT-5d52c7c8_S13_R1_001.postAdapt.postQual.fastq
python3 $pipeline_path -r $ref_18s -g $ref_gen -s $fastq