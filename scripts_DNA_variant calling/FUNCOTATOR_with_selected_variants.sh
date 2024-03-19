#!/bin/sh
for n in $(ls /mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/Variant_Calling/Select_variants/Lung33_tumor_S29.*.vcf.gz)
do
  echo $n
  sbatch --export=sample=$n,ALL FUNCOTATOR_with_selected_variants.sbatch
done
