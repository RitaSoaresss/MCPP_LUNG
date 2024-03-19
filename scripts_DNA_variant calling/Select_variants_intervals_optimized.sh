#!/bin/sh
for n in $(ls /mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/Variant_Calling/FilterMutectCalls/Lung33_tumor_S29.filtered.vcf.gz)
do
  echo $n
  sbatch --export=sample=$n,ALL Select_variants_intervals_optimized.sbatch
done
