#!/bin/bash
#SBATCH --job-name=2_MULTIQC_RNA
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --mem=240G
#SBATCH --cpus-per-task=40
#SBATCH --chdir=/home/rita.soares/scratch

path="/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/RNA/fastqc_rna"
multiqc_output="/mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/RNA/MULTIQC/multiqc_report_test.html"

 echo "START PROCESSING FASTQC"
 srun /mnt/beegfs/apptainer/images/ewels-multiqc-1.10.1.img $path -o $multiqc_output

#JOB STATISTICS (LIKE ELAPSED TIME AND MEMORY USAGE)
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS" -j $SLURM_JOB_ID
