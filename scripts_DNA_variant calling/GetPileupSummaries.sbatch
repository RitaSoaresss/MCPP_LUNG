#!/bin/bash
#SBATCH --job-name=GetPileupSummaries
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --mem=240G
#SBATCH --cpus-per-task=40
#SBATCH --chdir=/mnt/nfs/lobo/IMM-NFS/cfranco/scratch

variant="/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/gatk_files/small_exac_common_3.hg38.vcf.gz"
genomic_interval="/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/gatk_files/small_exac_common_3.hg38.vcf.gz"

# Run gatk MarkDuplicatesSpark
process_dir1="/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/ApplyBQSR_2"
process_dir2="/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/Variant_Calling/GetPileupSummaries"

    ##determine file prefix
    s=${sample}

    pf=$(echo ${s%.bqsr.bam*})
    s2=${pf#*/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/ApplyBQSR_2/}
    echo "s2 is $s2"
    echo "pf is $pf"
    BASE="$process_dir2/$s2"

 echo "APPLYBQSR file BAM is $BAM"
 BAM="$pf.bqsr.bam"

 echo "GetPileupSummaries OUTPUT will be $OUTPUT as SAMPLE.pileups.table"
 OUTPUT="$BASE.pileups.table"

 echo "srun /mnt/beegfs/apptainer/images/gatk_latest.sif gatk GetPileupSummaries -I $BAM -V $variant -L $genomic_interval -O $OUTPUT --tmp-dir /mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/TMP"
 srun /mnt/beegfs/apptainer/images/gatk_latest.sif gatk GetPileupSummaries -I $BAM -V $variant -L $genomic_interval -O $OUTPUT --tmp-dir /mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/TMP

# Job statistics (like elapsed time and memory usage)
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS" -j $SLURM_JOB_ID
