#!/bin/sh
for n in $(ls /mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/RNASeq/fastq_RNA/fastq_RNA/*_R2_001.fastq.gz)
do
 echo $n
 sbatch --export=sample=$n,ALL 1_FASTQC_RNA.sbatch
done
