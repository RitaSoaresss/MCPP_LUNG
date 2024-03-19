#!/bin/sh
for n in $(ls /mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/Variant_Calling/Mutect2/Lung64_tumor_S21.f1r2.tar.gz)
do
  echo $n
  sbatch --export=sample=$n,ALL LearnReadOrientationModel.sbatch
done
