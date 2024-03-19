#!/bin/bash
#SBATCH --job-name=Select_variants
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --mem=240G
#SBATCH --cpus-per-task=40
#SBATCH --chdir=/mnt/nfs/lobo/IMM-NFS/cfranco/scratch

# Replace input.vcf with the path to your FilterMutectCalls VCF file
# Replace output_prefix with the desired prefix for the output files

genome_reference="/mnt/nfs/lobo/SALMEIDA-NFS/lcosta/MCPP-Lung/bwa_prefiles/GRCh38.primary_assembly.genome.fa"

process_dir1="/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/Variant_Calling/FilterMutectCalls"
process_dir2="/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/Variant_Calling/Select_variants/chr_variants_9LC"

##determine file prefix
s=${sample}
pf=$(echo ${s%.filtered.vcf.gz*})
s2=${pf#*/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/Variant_Calling/FilterMutectCalls/}
echo "s2 is $s2"
echo "pf is $pf"
BASE="$process_dir2/$s2"

echo "INPUT OF VARIANTS FILE will be $VARIANTS_FILE as SAMPLE.filtered.vcf.gz"
VARIANTS_FILE="$pf.filtered.vcf.gz"

# Get the list of chromosomes dynamically
chromosomes=($(awk -v chrom_prefix="chr" '$1 ~ chrom_prefix {print $1}' "${genome_reference}.fai"))

for chr in "${chromosomes[@]}"; do
    echo "Processing chromosome: $chr"

    # Use awk to get the length of the chromosome from the reference genome index
    chromosome_length=$(awk -v chrom="$chr" '$1 == chrom {print $2}' "${genome_reference}.fai")

    echo "$chromosome_length"

    # Define the number of intervals you want
    num_intervals=4

    # Calculate the interval size
    interval_size=$((chromosome_length / num_intervals))

    # Generate intervals
    intervals=()
    for ((i = 0; i < num_intervals; i++)); do
        start=$((i * interval_size + 1))
        end=$((start + interval_size - 1))
        intervals+=("$start-$end")

        echo "srun /mnt/beegfs/apptainer/images/gatk_latest.sif gatk SelectVariants -R $genome_reference -V $VARIANTS_FILE -O $BASE.$chr.$start-$end.vcf.gz -L $chr:$start-$end"
        srun /mnt/beegfs/apptainer/images/gatk_latest.sif gatk SelectVariants -R $genome_reference -V $VARIANTS_FILE -O $BASE.$chr.$start-$end.vcf.gz -L $chr:$start-$end
    done
done

# Job statistics (like elapsed time and memory usage)
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS" -j $SLURM_JOB_ID

