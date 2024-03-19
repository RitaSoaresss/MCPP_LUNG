#!/bin/sh
for n in $(ls /mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/Variant_Calling/GetPileupSummaries/Lung64_tumor_S21.pileups.table)
do
  echo $n
  sbatch --export=sample=$n,ALL CalculateContamination.sbatch
done
