#!/bin/bash
#SBATCH --job-name=MergeBamAlignment_new_alignment
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --mem=240G
#SBATCH --cpus-per-task=40
#SBATCH --chdir=/home/rita.soares/scratch/

path="/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/NEW_ILLUMINA_PIPELINE"

# Define the directories for input unmapped bam, input unsorted bam and output merged bam
dir="$path/new_FastqToSam"
dir2="/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/NEW_ILLUMINA_DNA/new_samToBam"
dir3="$path/new_MergeBamAlignment"


# Merge unmapped bam and unsorted-aligned bam

echo "Start merging for new sample"

genome="/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/gatk/GRCh38.primary_assembly.genome.fa"


    ##determine file prefix
    s=${sample}
    pf3=$(echo ${s%.unmapped.*})
    echo "pf3 is $pf3"
    NAME1=$s
    BASE=${pf3#*/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/NEW_ILLUMINA_PIPELINE/new_FastqToSam/}
    NAME2=$dir2/$BASE.unsorted.bam

   echo "NAME1 is $NAME1"
   echo "NAME2 is $NAME2"
   echo "BASE is $BASE"

 echo "Start processing $NAME1 and $NAME2"

 echo "Getting merged bam..."
 echo "srun shifter java -XX:ParallelGCThreads=40 -jar /usr/picard/picard.jar MergeBamAlignment ALIGNED=$NAME2 UNMAPPED=$NAME1 O=$dir3/$BASE.merged.bam R=$genome TMP_DIR=/mnt/nfs/lobo/NMORAIS-NFS/imm-nas-p1A>"
 srun /mnt/beegfs/apptainer/images/picard_latest.sif java -XX:ParallelGCThreads=40 -jar /usr/picard/picard.jar MergeBamAlignment ALIGNED=$NAME2 UNMAPPED=$NAME1 O=$dir3/$BASE.merged.bam R=$genome SORT_ORDER=c>
# Job statistics (like elapsed time and memory usage)
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS" -j $SLURM_JOB_ID
