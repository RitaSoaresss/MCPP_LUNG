#!/bin/bash
#SBATCH --job-name=FilterMutectCalls
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --mem=240G
#SBATCH --cpus-per-task=40
#SBATCH --chdir=/mnt/nfs/lobo/IMM-NFS/cfranco/scratch

genome_reference="/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/bwa_prefiles/GRCh38.primary_assembly.genome.fa"


# Run gatk MarkDuplicatesSpark
process_dir1="/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/Variant_Calling/Mutect2"
process_dir2="/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/Variant_Calling/FilterMutectCalls"
process_dir3="/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/Variant_Calling/CalculateContamination"
process_dir4="/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/Variant_Calling/LearnReadOrientationModel"

    ##determine file prefix
    s=${sample}

    pf=$(echo ${s%.vcf.gz*})
    s2=${pf#*/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/Variant_Calling/Mutect2/}
    echo "s2 is $s2"
    echo "pf is $pf"
    BASE="$process_dir2/$s2"
    BASE_CONTAMINATION="$process_dir3/$s2"
    BASE_ARTIFACTS="$process_dir4/$s2"

 echo "INPUT OF VCF FILE CONTAINING VARIANTS WHICH IS THE OUTPUT FROM THE MUTECT2 STEP will be $VARIANT_FILE as SAMPLE.vcf.gz"
 VARIANT_FILE="$pf.vcf.gz"

 echo "CONTAMINATION_TABLE WHICH IS THE OUTPUT FROM THE CALCULATECONTAMINATION STEP will be $CONTAMINATION_TABLE as SAMPLE.contamination.table"
 CONTAMINATION_TABLE="$BASE_CONTAMINATION.contamination.table"

 echo "TABLES CONTAINING THE TUMOR SEGMENTS MINOR ALLELE FRACTIONS FOR GERMLINE HETS EMITTED BY THE CALCULATECONTAMINATION STEP will be $TUMOR_SEG as SAMPLE.tumor_segments.table"
 TUMOR_SEG="$BASE_CONTAMINATION.tumor_segments.table"

 echo "TABLES OF PRIOR ARTIFACT PROBABILITIES FOR THE READ ORIENTATION FILTER MODEL, ONE TABLE PER TUMOR SAMPLES will be $ob_priors as SAMPLE.artifact-prior.tar.gz"
 OB_PRIORS="$BASE_ARTIFACTS.artifact-prior.tar.gz"

 echo "THE MUTECT STATS FILE OUTPUT BY MUTECT2 will be $STATS_FILE as SAMPLE.vcf.gz.stats"
 STATS_FILE="$pf.vcf.gz.stats"

 echo "THE OUTPUT FILTERED VCF FILE will be $OUTPUT_FILTERED_VCF as SAMPLE.filtered.vcf.gz"
 OUTPUT_FILTERED_VCF="$BASE.filtered.vcf.gz"

 echo "srun /mnt/beegfs/apptainer/images/gatk_latest.sif gatk FilterMutectCalls -R $genome_reference -V $VARIANT_FILE --contamination-table $CONTAMINATION_TABLE --tumor-segmentation $TUMOR_SEG -O $OUTPUT_LEARN --tmp-dir /mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/TMP"
 srun /mnt/beegfs/apptainer/images/gatk_latest.sif gatk FilterMutectCalls -R $genome_reference -V $VARIANT_FILE --contamination-table $CONTAMINATION_TABLE --orientation-bias-artifact-priors $OB_PRIORS --stats $STATS_FILE --tumor-segmentation $TUMOR_SEG -O $OUTPUT_FILTERED_VCF --tmp-dir /mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/TMP

# Job statistics (like elapsed time and memory usage)
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS" -j $SLURM_JOB_ID
