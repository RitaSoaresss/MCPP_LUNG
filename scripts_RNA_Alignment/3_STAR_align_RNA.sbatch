#!/bin/bash
#SBATCH --job-name=3_STAR_align_RNA
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --mem=240G
#SBATCH --cpus-per-task=40
#SBATCH --chdir=/home/rita.soares/scratch
##determine file prefix

s=$sample
echo $s

t=(${s//"/"/ }) #criar lista
g=${s[${#s[@]}-1]} #selecionar elementos dessa lista
r=(${g//"_"/ }) #separa-los por "_"
l=${r[0]} #selcionar o 1º elemento
mkdir /mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/RNA/STAR_align/$l

pf=$(echo ${s%_R1_001*})
NAME1=$(echo ${pf}_R1_001.fastq.gz)
NAME2=$(echo ${pf}_R2_001.fastq.gz)
echo $pf
echo $NAME1
echo $NAME2

echo "Start processing $pf"

srun /mnt/beegfs/apptainer/images/star_latest.sif --genomeDir /mnt/nfs/lobo/SALMEIDA-NFS/lcosta/installations/annotations/star_index \
 --readFilesIn $NAME1 $NAME2 \
 --readFilesCommand gunzip -c \
 --runThreadN 40 \
 --outFilterMultimapScoreRange 1 \
 --outFilterMultimapNmax 10 \
 --outFilterMismatchNmax 10 \
 --alignIntronMax 500000 \
 --alignMatesGapMax 1000000 \
 --sjdbScore 2 \
 --alignSJDBoverhangMin 1 \
 --genomeLoad NoSharedMemory \
 --limitBAMsortRAM 70000000000 \
 --outFilterMatchNminOverLread 0.33 \
 --outFilterScoreMinOverLread 0.33 \
 --sjdbOverhang 100 \
 --outSAMstrandField intronMotif \
 --outSAMattributes NH HI NM MD AS XS \
 --sjdbGTFfile /mnt/nfs/lobo/SALMEIDA-NFS/lcosta/installations/annotations/gencode.v40lift37.annotation.gtf \
 --limitSjdbInsertNsj 1000000 \
 --outSAMunmapped None \
 --twopassMode Basic \
 --quantMode GeneCounts \
 --outFileNamePrefix /mnt/nfs/lobo/IMM-NFS/cfranco/MCPP-Lung/RNA/STAR_align/$l/"$l"_ \
 --outSAMmultNmax 1 \

echo "... done processing $pf"

#JOB STATISTICS (LIKE ELAPSED TIME AND MEMORY USAGE)
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS" -j $SLURM_JOB_ID
