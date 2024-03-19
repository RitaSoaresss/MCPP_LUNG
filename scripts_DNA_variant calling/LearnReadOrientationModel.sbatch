#!/bin/bash
#SBATCH --job-name=LearnReadOrientationModel
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --mem=240G
#SBATCH --cpus-per-task=40
#SBATCH --chdir=/mnt/nfs/lobo/IMM-NFS/cfranco/scratch

# Run gatk MarkDuplicatesSpark
process_dir1="/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/Variant_Calling/Mutect2"
process_dir2="/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/Variant_Calling/LearnReadOrientationModel"

    ##determine file prefix
    s=${sample}

    pf=$(echo ${s%.f1r2.tar.gz*})
    s2=${pf#*/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/Variant_Calling/Mutect2/}
    echo "s2 is $s2"
    echo "pf is $pf"
    BASE="$process_dir2/$s2"

 echo "INPUT OF THE F1R2 COUNTS will be $F1R2_COUNTS as SAMPLE.f1r2.tar.gz"
 F1R2_COUNTS="$pf.f1r2.tar.gz"

 echo "OUTPUT will be $OUTPUT_LEARN as SAMPLE.artifact-prior.tar.gz"
 OUTPUT_LEARN="$BASE.artifact-prior.tar.gz"

 echo "srun /mnt/beegfs/apptainer/images/gatk_latest.sif gatk LearnReadOrientationModel -I $F1R2_COUNTS -O $OUTPUT_LEARN --tmp-dir /mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/TMP"
 srun /mnt/beegfs/apptainer/images/gatk_latest.sif gatk LearnReadOrientationModel -I $F1R2_COUNTS -O $OUTPUT_LEARN --tmp-dir /mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/TMP

# Job statistics (like elapsed time and memory usage)
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS" -j $SLURM_JOB_ID
