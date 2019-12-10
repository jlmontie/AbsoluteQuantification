INDEX=GCF_000001405.39_GRCh38.p13_rna_from_genomic.idx
OUTDIR=IDBD-D100376
FQDIR=/home/jmontgomery/mnt/arup_urine_fqs
FQ=190828-1-1-IDBD-D100376-d-01-AHMW75AFXY-TGATTATACG-TGTAATCGAC-5d6e759f_S1_R1_001.postAdapt.postQual.fastq
docker run -v $PWD:/mnt -v $FQDIR:/mnt2 kallisto-img kallisto quant --single -l 300 -s 30 -i /mnt/$INDEX -o /mnt/$OUTDIR /mnt2/$FQ