#!/bin/sh
for n in $(ls /mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/NEW_ILLUMINA_DNA/new_samToBam/*.unsorted.bam)
do
  echo $n
  sbatch --export=sample=$n,ALL SortSamSpark_new_illumina.sbatch
done