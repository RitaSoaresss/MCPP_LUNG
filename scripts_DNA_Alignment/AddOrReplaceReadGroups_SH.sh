#!/bin/sh
for n in $(ls /mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/NEW_ILLUMINA_PIPELINE/new_MergeBamAlignment/*.merged.bam)
do
  echo $n
  sbatch --export=sample=$n,ALL AddOrReplaceReadGroups_new_illumina.sbatch
done