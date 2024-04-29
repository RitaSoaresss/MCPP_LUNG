#!/bin/bash
#SBATCH --job-name=MULTIQC_DNA
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --mem=240G
#SBATCH --cpus-per-task=40
#SBATCH --chdir=/home/rita.soares/scratch

path="/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/fastqc_DNA"
multiqc_output="/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/DNA/MULTIQC_DNA/multiqc_report_dna.html"

 echo "START PROCESSING FASTQC"
 srun /mnt/beegfs/apptainer/images/ewels-multiqc-1.10.1.img $path -o $multiqc_output

#JOB STATISTICS (LIKE ELAPSED TIME AND MEMORY USAGE)
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS" -j $SLURM_JOB_ID
