#!/bin/sh
for n in $(ls /mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/NEW_ILLUMINA_PIPELINE/new_fastp/*_R1.trimmed.fastq.gz)
do
  echo $n
  sbatch --export=sample=$n,ALL bwa_new_illumina.sbatch
done