#!/bin/sh
for n in $(ls /mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/SortSam_MarkDuplicates/*.markeduplicates.sorted.bam)
do
  echo $n
  sbatch --export=sample=$n,ALL BaseRecalibrator.sbatch
done
