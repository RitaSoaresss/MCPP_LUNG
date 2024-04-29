#!/bin/bash
#SBATCH --job-name=bwa.Group1
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --mem=240G
#SBATCH --cpus-per-task=40
#SBATCH --chdir=/home/rita.soares/scratch

path="/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/NEW_ILLUMINA_PIPELINE"
genome2="/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/bwa_prefiles/GRCh38.primary_assembly.genome.fa"

# Align trimmed fastq files to the GRCh38 genome
process_dir="$path/new_fastp"
process_dir2="$path/new_BWA"

    ##determine file prefix
    s=${sample}

    pf=$(echo ${s%_R1.trimmed*})
    NAME1=$(echo ${pf}_R1.trimmed.fastq.gz)
    NAME2=$(echo ${pf}_R2.trimmed.fastq.gz)

    s2=${pf#*/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/NEW_ILLUMINA_PIPELINE/new_fastp}
    echo "s2 is $s2"
    echo "pf is $pf"
    BASE="$process_dir2/$s2"

    echo "BASE file name is $BASE"
    echo "NAME1 IS $NAME1"
    echo "NAME2 IS $NAME2"

 echo "Start processing $pf"

 echo "Aligning to genome..."
 echo "bwa mem -t 40 $genome2 $NAME1 $NAME2 > $BASE.sam"
 srun /mnt/beegfs/apptainer/images/bwa_latest.sif bwa mem -t 40 $genome2 $NAME1 $NAME2 > $BASE.sam

# Job statistics (like elapsed time and memory usage)
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS" -j $SLURM_JOB_ID