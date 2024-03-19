#!/bin/sh
for n in $(ls /mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/Variant_Calling/Mutect2/Lung64_tumor_S21.vcf.gz)
do
  echo $n
  sbatch --export=sample=$n,ALL FilterMutectCalls.sbatch
done
