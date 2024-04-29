#!/bin/sh
for n in $(ls /mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/NEW_ILLUMINA_PIPELINE/new_bwa/*.sam)
do
  echo $n
  sbatch --export=sample=$n,ALL samToBam-Group1.sbatch
done