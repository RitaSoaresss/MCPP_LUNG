#!/bin/sh
for n in $(ls /mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/ApplyBQSR_2/Lung64_tumor_S21.bqsr.bam)
do
  echo $n
  sbatch --export=sample=$n,ALL GetPileupSummaries.sbatch
done
