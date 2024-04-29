#!/bin/bash
#SBATCH --job-name=bwa.Group1
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --mem=240G
#SBATCH --cpus-per-task=40
#SBATCH --image=biocontainers/samtools:v1.9-4-deb_cv1
#SBATCH --chdir=/home/rita.soares/scratch/
#SBATCH --exclude=compute-[1,5,9,13,17,21]

path="/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung"

# Transform bwa-produced .sam file to .bam format
process_dir="$path/bwa"
process_dir2="$path/samToBam"

    ##determine file prefix
    s=${sample}

    pf=$(echo ${s%.sam*})
    s2=${pf#*/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/bwa/}
    echo "s2 is $s2"
    echo "pf is $pf"
    BASE="$process_dir2/$s2"

 echo "convert sam to bam"
 SAM="$pf.sam"

 echo "SAM is $SAM"
 echo "samtools view -S -b $SAM > "$BASE.unsorted.bam""
 srun shifter samtools view -S -b $SAM > "$BASE.unsorted.bam"

 #echo "Done bam file"

# Job statistics (like elapsed time and memory usage)
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS" -j $SLURM_JOB_ID

