#!/bin/bash
for f in $(/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/RNA/fastq_RNA/Lung50_RNA_tumorRNA_S6_R1*)
do
 echo $f
 sbatch --export=sample=$f,ALL STAR_2.7.9a_2pass_pairend_end.trim.sbatch
done
