#!/bin/sh
for n in $(ls /mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/AddOrReplaceReadGroups/Lung*)
do
  echo $n
  sbatch --export=sample=$n,ALL MarkDuplicates.sbatch
done