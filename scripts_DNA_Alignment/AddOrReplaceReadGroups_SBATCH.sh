#!/bin/bash
#SBATCH --job-name=AddOrReplaceReadGroups_new_illumina
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --mem=240G
#SBATCH --cpus-per-task=40
#SBATCH --chdir=/home/rita.soares/scratch/

path="/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/NEW_ILLUMINA_PIPELINE"

# Run gatk SortSamSpark
process_dir="$path/new_MergeBamAlignment"
process_dir2="$path/new_AddOrReplaceRead"

    ##determine file prefix
    s=${sample}

    pf=$(echo ${s%.merged.bam*})
    s2=${pf#*/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/NEW_ILLUMINA_PIPELINE/new_MergeBamAlignment/}
    echo "s2 is $s2"
    echo "pf is $pf"
    BASE="$process_dir2/$s2"

 echo "AddOrReplaceReadGroups (Picard)"
 READGROUPS="$BASE.readgroups.bam"

 echo "ReadGroups bam will be $READGROUPS"
 echo "srun shifter java -XX:ParallelGCThreads=40 -jar /usr/picard/picard.jar AddOrReplaceReadGroups I=$s O=$READGROUPS RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=2"
 srun /mnt/beegfs/apptainer/images/picard_latest.sif java -XX:ParallelGCThreads=40 -jar /usr/picard/picard.jar AddOrReplaceReadGroups I=$s O=$READGROUPS RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=2


# Job statistics (like elapsed time and memory usage)
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS" -j $SLURM_JOB_ID