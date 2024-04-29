#!/bin/bash
#SBATCH --job-name=FASTQC_RNA_2
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --mem=240G
#SBATCH --cpus-per-task=40
#SBATCH --chdir=/home/rita.soares/scratch

# MOVE THE RESULT FROM FASTQC TO fastqc_rna

process_dir="/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/RNASeq/fastq_RNA/fastq_RNA"
process_dir2="/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/RNA/fastqc_rna"

        ##DETERMINE THE FILE PREFIX
        s=${sample}

#       pf=$(echo ${s%_R1_001.*})
#       NAME1=$(echo ${pf}_R1_001.fastq.gz)

 echo "START PROCESSING FASTQC"
 srun /mnt/beegfs/apptainer/images/biocontainers-fastqc-v0.11.9_cv8.img fastqc $s -o $process_dir2

#JOB STATISTICS (LIKE ELAPSED TIME AND MEMORY USAGE)
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS" -j $SLURM_JOB_ID
