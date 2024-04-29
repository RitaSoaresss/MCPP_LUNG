#!/bin/bash
#SBATCH --job-name=ApplyBQSR
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --mem=240G
#SBATCH --cpus-per-task=40
#SBATCH --chdir=/mnt/nfs/lobo/IMM-NFS/cfranco/scratch

path="/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung"

genome_reference="/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/bwa_prefiles/GRCh38.primary_assembly.genome.fa"

# Run gatk MarkDuplicatesSpark
process_dir1="/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/SortSam_MarkDuplicates"
process_dir2="/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/BaseRecalibrator_snps_indels"
process_dir3="/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/ApplyBQSR"

    ##determine file prefix
    s=${sample}

    pf=$(echo ${s%.markeduplicates.sorted.bam*})
    s2=${pf#*/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/SortSam_MarkDuplicates/}
    echo "s2 is $s2"
    echo "pf is $pf"
    BASE="$process_dir2/$s2"
    BASE_OUTPUT="$process_dir3/$s2"

 echo "Mark Duplicates BAM is $BAM"
 BAM="$pf.markeduplicates.sorted.bam"

 echo "recal.table from the BaseRecalibrator is $bqsr_recal_file"
 bqsr_recal_file="$BASE.recal.table"

 echo "ApplyBQSR OUTPUT will be $OUTPUT as SAMPLE_bqsr.bam"
 OUTPUT="$BASE_OUTPUT.bqsr.bam"

 echo "srun /mnt/beegfs/apptainer/images/gatk_latest.sif gatk ApplyBQSR -R $genome_reference -I $BAM --bqsr-recal-file $bqsr_recal_file -O $OUTPUT --tmp-dir /mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/TMP"
 srun /mnt/beegfs/apptainer/images/gatk_latest.sif gatk ApplyBQSR -R $genome_reference -I $BAM --bqsr-recal-file $bqsr_recal_file -O $OUTPUT --tmp-dir /mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/TMP

# Job statistics (like elapsed time and memory usage)
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS" -j $SLURM_JOB_ID
