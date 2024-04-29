#!/bin/sh
for n in $(ls /mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/SortSam_MarkDuplicates/Lung54_tumor_S12.markeduplicates.sorted.bam)
do
  echo $n
  sbatch --export=sample=$n,ALL ApplyBQSR.sbatch
done
