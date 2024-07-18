# RNA & DNA PRE-PROCESSING

## A. RNA PRE-PROCESSING AND DECONVOLUTION

### A.1 RNA PRE-PROCESSING WORKFLOW:

![image](https://github.com/user-attachments/assets/2f9b1152-8812-4f57-975b-39b6333ff15d)

Here's a list of the RNA pre-processing scripts and a brief description of them in order:

#### 1_FASTQC_RNA

The FastQ files provided from the whole sequencing, pass first through a quality control analysis using FASTQC tool

#### 2_MULTIQC_RNA

The MULTIQC tool consists in providing a html report of all the samples fastqc files into one.

#### 3_STAR_align_RNA

Aligns the clean fastq files to the reference genome GRCh38 (depends on the type of data), to output a sam format file, as well as other files, such as, “ReadsPerGene.out.tab" used then in the differential gene expression with DESEQ2; ESTIMATE analysis ; and deconvolution analysis with CIBERSORTx.

### A.2 DECONVOLUTION WORKFLOW:

![image](https://github.com/user-attachments/assets/15ec7b08-3f66-441a-b5c6-d84b97505925)

CIBERSORTx is a deconvolution web tool (https://cibersortx.stanford.edu/), that was performed with a single cell public dataset from the single cell atlas (https://www.ebi.ac.uk/gxa/sc/home). This dataset was used as reference purified gene expression profile of the different cell types, to origine a signature/barcode matrix. Single cell public dataset E-MTAB-6653 (https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-6653) from three patients with untreated NSCLC as signature matrix.

#### SC_LUNG_3-lung-carcinomas.Rmd

The annotation of the single cell data was done by the Seurat package in R (https://satijalab.org/seurat/). The workflow of the annotation divides into 5 major steps: Quality control and subset of single cell data; Normalization of the data; Identification of highly variable features (feature selection) and Scaling of the data; Linear dimensional reduction and clustering; Cell annotation and diagnosis.

#### CIBERSORTx

## B. DNA PRE-PROCESSING WORKFLOW:

The DNA pre-processing is divided into 2 parts: 

### B.1 - Variant Alignment

  ![image](https://github.com/user-attachments/assets/c107d7c3-0091-4361-a8ba-b6fbd9ce882e)

  Here's a list of the DNA Variant Alignment scripts and a brief description of them in order:

  #### 1_FASTP and 1_MULTIQC

  #### 2_BWA_MEM

  #### 3_SAMTools

  #### 4_FastqToSam

  #### 5_MergeBamAlinment

  #### 6_MarkDuplicates

  #### 7_SortSam

  #### 8_BaseRecalibrator

  #### 9_ApplyBQSR

### B.2 - Variant Calling

  ![image](https://github.com/user-attachments/assets/a3c5191b-b77b-43e9-b6c5-3c4c4347f479)

  Here's a list of the DNA Variant Calling scripts and a brief description of them in order:

  #### 1_Mutect2

  #### 2_GetPileupSummaries

  #### 3_CalculateContamination

  #### 4_LearnReadOrientationModel

  #### 5_FilterMutectCalls

  #### 6_SelectVariants

  #### 7_Funcotator

