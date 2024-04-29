#!/bin/bash
#SBATCH --job-name=BaseRecalibrator
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --mem=240G
#SBATCH --cpus-per-task=40
#SBATCH --chdir=/mnt/nfs/lobo/IMM-NFS/cfranco/scratch

path="/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung"

genome_reference="/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/gatk/GRCh38.primary_assembly.genome.fa"
known_sites="/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/gatk_files/Homo_sapiens_assembly38.dbsnp138.vcf"
known_sites_indels="/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/gatk_files/Homo_sapiens_assembly38.known_indels.vcf.gz"

# Run gatk MarkDuplicatesSpark
process_dir="/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/SortSam_MarkDuplicates"
process_dir2="/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/BaseRecalibrator_snps_indels"

    ##determine file prefix
    s=${sample}

    pf=$(echo ${s%.markeduplicates.sorted.bam*})
    s2=${pf#*/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/SortSam_MarkDuplicates/}
    echo "s2 is $s2"
    echo "pf is $pf"
    BASE="$process_dir2/$s2"

 echo "Mark Duplicates BAM is $BAM"
 BAM="$pf.markeduplicates.sorted.bam"

 echo "BaseRecalibrator OUTPUT will be $OUTPUT as SAMPLE.recal.table"
 OUTPUT="$BASE.recal.table"

 echo "srun /mnt/beegfs/apptainer/images/gatk_latest.sif  gatk BaseRecalibrator -I $BAM -R $genome_reference --known-sites $known_sites --known-sites $known_sites_indels -O $OUTPUT --TMP_DIR /mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/TMP"
 srun /mnt/beegfs/apptainer/images/gatk_latest.sif  gatk BaseRecalibrator -I $BAM -R $genome_reference --known-sites $known_sites --known-sites $known_sites_indels -O $OUTPUT --tmp-dir /mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/TMP

# Job statistics (like elapsed time and memory usage)
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS" -j $SLURM_JOB_ID
