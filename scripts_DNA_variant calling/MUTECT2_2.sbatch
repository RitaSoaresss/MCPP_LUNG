#!/bin/bash
#SBATCH --job-name=MUTECT2_2
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --mem=240G
#SBATCH --cpus-per-task=40
#SBATCH --chdir=/mnt/nfs/lobo/IMM-NFS/cfranco/scratch

genome_reference="/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/bwa_prefiles/GRCh38.primary_assembly.genome.fa"
germline_file="/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/gatk/af-only-gnomad.hg38.vcf"
painel_of_normals="/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/gatk/1000g_pon.hg38.vcf"

# Run gatk MarkDuplicatesSpark
process_dir1="/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/ApplyBQSR_2"
process_dir2="/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/Variant_Calling/Mutect2"

    ##determine file prefix
    s=${sample}

    pf=$(echo ${s%.bqsr.bam*})
    s2=${pf#*/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/ApplyBQSR_2/}
    echo "s2 is $s2"
    echo "pf is $pf"
    BASE="$process_dir2/$s2"

 echo "APPLYBQSR file BAM is $BAM"
 BAM="$pf.bqsr.bam"

 echo "F1R2 counts output file is $F1R2_COUNTS"
 F1R2_COUNTS="$BASE.f1r2.tar.gz"

 echo "MUTECT2 OUTPUT will be $OUTPUT as SAMPLE.vcf.gz"
 OUTPUT="$BASE.vcf.gz"

 echo "srun /mnt/beegfs/apptainer/images/gatk_latest.sif gatk Mutect2 -R $genome_reference -I $BAM --f1r2-tar-gz $F1R2_COUNTS --germline-resource $germline_file --panel-of-normals $painel_of_normals -O $OUTPUT --tmp-dir /mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/TMP"
 srun /mnt/beegfs/apptainer/images/gatk_latest.sif gatk Mutect2 -R $genome_reference -I $BAM --f1r2-tar-gz $F1R2_COUNTS --germline-resource $germline_file --panel-of-normals $painel_of_normals -O $OUTPUT --tmp-dir /mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/TMP

# Job statistics (like elapsed time and memory usage)
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS" -j $SLURM_JOB_ID
