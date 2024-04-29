#!/bin/sh
for n in $(ls /mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Sarcoma/raw_data/New_Illumina_sent/Lung/*_R1_001.fastq.gz)
do
  echo $n
  sbatch --export=sample=$n,ALL fastp_new_illumina.sbatch
done