#!/bin/bash
#SBATCH --job-name=SortSamSpark_new_illumina
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --mem=240G
#SBATCH --cpus-per-task=40
#SBATCH --chdir=/home/rita.soares/scratch/

path="/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/NEW_ILLUMINA_PIPELINE"

# Run gatk SortSamSpark
process_dir="/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/NEW_ILLUMINA_DNA/new_samToBam"
process_dir2="$path/new_sortedbam"

    ##determine file prefix
    s=${sample}

    pf=$(echo ${s%.unsorted.bam*})
    s2=${pf#*/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/NEW_ILLUMINA_DNA/new_samToBam/}
    echo "s2 is $s2"
    echo "pf is $pf"
    BASE="$process_dir2/$s2"

 echo "sort bam"
 BAM="$pf.unsorted.bam"
 SORTED="$BASE.sorted.bam"

 echo "BAM is $BAM"
 echo "Sorted BAM will be $SORTED"

 echo "srun /mnt/beegfs/apptainer/images/gatk_latest.sif gatk SortSamSpark -I $BAM -O $SORTED --  --spark-runner LOCAL --spark-master 'local[*]' --tmp-dir /mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/tmp/"
 srun /mnt/beegfs/apptainer/images/gatk_latest.sif gatk SortSamSpark -I $BAM -O $SORTED --  --spark-runner LOCAL --spark-master 'local[*]' --tmp-dir /mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/tmp/
# Job statistics (like elapsed time and memory usage)
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS" -j $SLURM_JOB_ID