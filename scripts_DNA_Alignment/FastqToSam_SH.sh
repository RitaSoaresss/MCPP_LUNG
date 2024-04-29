#!/bin/sh
for n in $(ls /mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/NEW_ILLUMINA_PIPELINE/new_fastp/10LC_tumor_S2_R1.trimmed.fastq.gz)
do
  echo $n
  sbatch --export=sample=$n,ALL FastqToSam_10_new_illumina.sbatch
done