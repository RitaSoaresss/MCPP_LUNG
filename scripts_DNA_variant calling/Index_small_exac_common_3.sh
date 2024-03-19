#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --mem=240G
#SBATCH --cpus-per-task=40
#SBATCH --chdir=/mnt/nfs/lobo/IMM-NFS/cfranco/scratch

INPUT="/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/gatk_files/small_exac_common_3.hg38.vcf.gz"

srun /mnt/beegfs/apptainer/images/gatk_latest.sif gatk IndexFeatureFile -I $INPUT

# Job statistics (like elapsed time and memory usage)
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS" -j $SLURM_JOB_ID
