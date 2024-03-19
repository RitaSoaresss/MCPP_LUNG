#!/bin/bash
#SBATCH --job-name=FUNCOTATOR_with_selected_variants
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --mem=240G
#SBATCH --cpus-per-task=40
#SBATCH --chdir=/mnt/nfs/lobo/IMM-NFS/cfranco/scratch

genome_reference="/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/gatk_files/GRCh38.primary_assembly.genome.fa"
DATASOURCE="/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/gatk_files/funcotator_dataSources.v1.7.20200521s"

# Run gatk MarkDuplicatesSpark
process_dir1="/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/Variant_Calling/Select_variants"
process_dir2="/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/Variant_Calling/Select_variants/FUNCOTATOR_selected_variants_33LC"

    ##determine file prefix
    s=${sample}

    pf=$(echo ${s%.vcf.gz*})
    s2=${pf#*/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/Variant_Calling/Select_variants/}
    echo "s2 is $s2"
    echo "pf is $pf"
    BASE="$process_dir2/$s2"

 echo "INPUT OF VARIANTS FILE will be $VARIANTS_FILE as SAMPLE.filtered.vcf.gz"
 VARIANTS_FILE="$pf.vcf.gz"

 echo "OUTPUT will be $OUTPUT_LEARN as SAMPLE.FiltPASS.recode_funcotator.vcf"
 OUTPUT_MAF="$BASE.FiltPASS.maf"

 echo "srun /mnt/beegfs/apptainer/images/gatk_latest.sif gatk Funcotator -R $genome_reference -V $VARIANTS_FILE -O $OUTPUT_MAF --output-file-format MAF --ref-version hg38 --data-sources-path $DATASOURCE --disable-sequence-dictionary-validation"
 srun /mnt/beegfs/apptainer/images/gatk_latest.sif gatk Funcotator -R $genome_reference -V $VARIANTS_FILE -O $OUTPUT_MAF --output-file-format MAF --ref-version hg38 --data-sources-path $DATASOURCE --disable-sequence-dictionary-validation

# Job statistics (like elapsed time and memory usage)
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS" -j $SLURM_JOB_ID
