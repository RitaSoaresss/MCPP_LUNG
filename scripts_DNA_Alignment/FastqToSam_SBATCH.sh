#!/bin/bash
#SBATCH --job-name=FastqToSam_10_new_illumina
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --mem=240G
#SBATCH --cpus-per-task=40
#SBATCH --chdir=/home/rita.soares/scratch/

path="/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/NEW_ILLUMINA_PIPELINE"

# Run picard FastqToSam
process_dir="$path/new_fastp"
process_dir2="$path/new_FastqToSam"

##determine file prefix
    s=${sample}

    pf=$(echo ${s%_R1.trimmed*})
    NAME1=$(echo ${pf}_R1.trimmed.fastq.gz)
    NAME2=$(echo ${pf}_R2.trimmed.fastq.gz)

    s2=${pf#*/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/new_fastp/}
    echo "s2 is $s2"
    echo "pf is $pf"
    BASE="$process_dir2"

    echo "BASE file name is $BASE"
    echo "NAME1 IS $NAME1"
    echo "NAME2 IS $NAME2"

 echo "Start processing $pf"
 echo "FastqToSam (Picard) producing uBAM"
 uBAM="$BASE/10LC_tumor_S2.unmapped.bam"

 echo "uBAM will be $uBAM"
 echo "srun shifter java -XX:ParallelGCThreads=40 -jar /usr/picard/picard.jar FastqToSam F1=$NAME1 F2=$NAME2 O=$uBAM SM=$s2 TMP_DIR=/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/tmp/"
 srun /mnt/beegfs/apptainer/images/picard_latest.sif java -XX:ParallelGCThreads=40 -jar /usr/picard/picard.jar FastqToSam F1=$NAME1 F2=$NAME2 O=$uBAM SM=$s2 TMP_DIR=/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lun>

# Job statistics (like elapsed time and memory usage)
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS" -j $SLURM_JOB_ID