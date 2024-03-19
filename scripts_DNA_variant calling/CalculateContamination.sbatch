#!/bin/bash
#SBATCH --job-name=CalculateContamination
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --mem=240G
#SBATCH --cpus-per-task=40
#SBATCH --chdir=/mnt/nfs/lobo/IMM-NFS/cfranco/scratch

# Run gatk MarkDuplicatesSpark
process_dir1="/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/Variant_Calling/GetPileupSummaries"
process_dir2="/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/Variant_Calling/CalculateContamination"

    ##determine file prefix
    s=${sample}

    pf=$(echo ${s%.pileups.table*})
    s2=${pf#*/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/Variant_Calling/GetPileupSummaries/}
    echo "s2 is $s2"
    echo "pf is $pf"
    BASE="$process_dir2/$s2"

 echo "GetPileupSummaries table file is $PILEUPS"
 PILEUPS="$pf.pileups.table"

 echo "OUTPUT_SEG will be $TUMOR_SEG_TABLE as SAMPLE.tumor_segments.table"
 TUMOR_SEG_TABLE="$pf.tumor_segments.table"

 echo "CALCULATE CONTAMINATION OUTPUT will be $OUTPUT_CONTAMINATION as SAMPLE.contamination.table"
 OUTPUT_CONTAMINATION="$BASE.contamination.table"

 echo "srun /mnt/beegfs/apptainer/images/gatk_latest.sif gatk CalculateContamination -I $PILEUPS --tumor-segmentation $TUMOR_SEG_TABLE -O $OUTPUT --tmp-dir /mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/TMP"
 srun /mnt/beegfs/apptainer/images/gatk_latest.sif gatk CalculateContamination -I $PILEUPS --tumor-segmentation $TUMOR_SEG_TABLE -O $OUTPUT_CONTAMINATION --tmp-dir /mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/TMP

# Job statistics (like elapsed time and memory usage)
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS" -j $SLURM_JOB_ID
