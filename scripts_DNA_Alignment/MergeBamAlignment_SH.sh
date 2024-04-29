#!/bin/sh
for n in $(ls /mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/NEW_ILLUMINA_PIPELINE/new_FastqToSam/26LC_tumor_S7.unmapped.bam)
do
  echo $n
  sbatch --export=sample=$n,ALL MergeBamAlignment_new_illumina.sbatch
done