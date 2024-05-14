rm(list = ls())

#GENE NAMES:
library(readr)
BULKDATA <- read_delim("C:/Users/ritas/Desktop/TESE_IMM/cibersort_my_lung_samples/reads_per_gene_star_MCPP-LUNG/BULKDATA.txt", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)
View(BULKDATA)
BULKDATA$Lung27<-NULL
BULK_DATA<- as.data.frame(BULKDATA[!duplicated(BULKDATA$Gene), ])
BULK_DATA_rownames <- data.frame(BULK_DATA[,-1], row.names=BULK_DATA[,1])
BULK_DATA_rownames[is.na(BULK_DATA_rownames)] <- 0
View(BULK_DATA_rownames)

#Create the sample sheet to define the group variable (1º/2º/3º treatment line):
library(readxl)
SAMPLES_GROUPS_TREATMENT_RNA <- read_excel("C:/Users/ritas/Desktop/TESE_IMM/cibersort_my_lung_samples/SAMPLES_GROUPS_TREATMENT_BIOPSIE_RNA.xlsx")
View(SAMPLES_GROUPS_TREATMENT_RNA)
SAMPLES_GROUPS_TREATMENT_RNA<-as.data.frame(SAMPLES_GROUPS_TREATMENT_RNA)

#transform the samples names into rownames
SAMPLES_GROUPS_TREATMENT_RNA_<- data.frame(SAMPLES_GROUPS_TREATMENT_RNA[,-1], row.names=SAMPLES_GROUPS_TREATMENT_RNA[,1])
colnames(SAMPLES_GROUPS_TREATMENT_RNA_)<-c("TREATMENT_GROUPS","TYPE","RESPONSE","TYPE_PROGRESSIVE")
View(SAMPLES_GROUPS_TREATMENT_RNA_)

#remove the patients that don't have any information about their type of biopsie and also no available information regarding response:
metadata<-SAMPLES_GROUPS_TREATMENT_RNA[complete.cases(SAMPLES_GROUPS_TREATMENT_RNA$TYPE),]
View(metadata)#34 de 44 patients

rownames(metadata)<-metadata$SAMPLES
#put the samples in the same order in the count data as in the col data:
BULK_DATA_rownames_<-BULK_DATA_rownames[,rownames(metadata)]

#1º CONFIRM IF THE ROW NAMES IN BULK_DATA_ROWNAMES (COUNTS_DATA) ARE ALL PRESENT IN SAMPLES_GROUPS_TREATMENT_RNA (COLDATA)
all(colnames(BULK_DATA_rownames_) %in%rownames(metadata))#TRUE

#2º CONFIRM IF THOSE SAMPLE NAMES ARE IN THE SAME ORDER IN BOTH DATAFRAMES:
all(colnames(BULK_DATA_rownames_)==rownames(metadata))

View(BULK_DATA_rownames_)
BULK_DATA_rownames_$Lung9<-NULL

metadata <- subset(metadata, rownames(metadata) != "Lung9")
View(metadata)

#############################################################################
#FILTER GENES- GET RIDE OF THE PSEUDOGENES AND NO-CODING GENES:
#############################################################################
View(BULK_DATA_rownames_)#39,764 genes in total

library(biomaRt)  
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'useast.ensembl.org')
genes <- biomaRt::getBM(attributes = c("external_gene_name", "chromosome_name","transcript_biotype"), filters = c("transcript_biotype","chromosome_name"),values = list("protein_coding",c(1:22,"M","X","Y")), mart = mart)
coding_genes_until_Y<-genes
View(coding_genes_until_Y)#19,366 GENES

#GET THE LIST OF THOSE GENES:
coding_genes_list<-coding_genes_until_Y$external_gene_name
View(coding_genes_list)

rownames(BULK_DATA_rownames_)

#ensemble part:

#A.CREATE THE DATAFRAME COLUMN WITH THE ENSEMBLE NAMES AND ADD TO THE ORIGINAL DATAFRAME:
d<-as.data.frame(rownames(BULK_DATA_rownames_))
d=sub("\\..*", "", rownames(BULK_DATA_rownames_))
d=as.data.frame(d)
View(d)
b<-cbind(d,BULK_DATA_rownames_)
View(b)

#B.SELECT THE Keep the first occurrence of each duplicated row
# Identify duplicated rows
duplicates <- duplicated(b$d)#20 duplicated

# Keep the first occurrence of each duplicated row
my_df <- b[!duplicates, ]

#CONFIRM IF THERE ARE ANY DUPLICATES IN THE FINAL DATAFRAME:
c<-duplicated(my_df)
table(c)#FALSE-40948 - SO ALL ARE UNIQUE

final_db<-my_df
rownames(final_db)<-final_db$d
final_db$d<-NULL
View(final_db)#39764 genes in total

BULK_DATA_rownames_<-final_db
##########################

BULK_DATA_rownames_df <- BULK_DATA_rownames_[row.names(BULK_DATA_rownames_) %in% coding_genes_list, ]
View(BULK_DATA_rownames_df)#19,187 genes

#put the samples in the same order in the count data as in the col data:
BULK_DATA_rownames_df_<-BULK_DATA_rownames_df[,rownames(metadata)]

#1º CONFIRM IF THE ROW NAMES IN BULK_DATA_ROWNAMES (COUNTS_DATA) ARE ALL PRESENT IN SAMPLES_GROUPS_TREATMENT_RNA (COLDATA)
all(colnames(BULK_DATA_rownames_df_) %in%rownames(metadata))#TRUE

#2º CONFIRM IF THOSE SAMPLE NAMES ARE IN THE SAME ORDER IN BOTH DATAFRAMES:
all(colnames(BULK_DATA_rownames_df_)==rownames(metadata))

##################################################################################################################

library(DESeq2)
STAR_LUNG_deseq2<- DESeqDataSetFromMatrix(countData=BULK_DATA_rownames_df_,
                                          colData = metadata,
                                          design=~TYPE)
STAR_LUNG_deseq2#43 samples primary vs metastases
View(STAR_LUNG_deseq2)

##############################FILTERING BEFORE FITTING MODEL###########3
# Number of genes before filtering:
nrow(STAR_LUNG_deseq2)#19187 genes with conversion/63368 with ensembl symbol

# Filter
STAR_LUNG_filtered<- STAR_LUNG_deseq2[rowSums(counts(STAR_LUNG_deseq2)) > 5, ]


# Number of genes left after low-count filtering:
nrow(STAR_LUNG_filtered)#18584 option C gene names/ 18800 genes with ensemnl symbol

library(DESeq2)
ddsHTSeq__LUNG_filtered <- DESeq(STAR_LUNG_filtered)
ddsHTSeq__LUNG_filtered
View(ddsHTSeq__LUNG_filtered)
names(ddsHTSeq__LUNG_filtered)

#extract the results in the form of a table using the results function:
res<- results(ddsHTSeq__LUNG_filtered)
View(res)
View(res)
head(res)
res<-na.omit(res)
res<-res[res$baseMean>50,]

summary(signficant_res)
##############################################################################

signficant_res<-subset(res,padj<0.05)
View(signficant_res)

signficant_res_dataframe<-as.data.frame(signficant_res)
signficant_res_with_rownames <- cbind(rownames(signficant_res_dataframe), signficant_res_dataframe)
colnames(signficant_res_with_rownames)[1] <- "GeneName"

write.xlsx(signficant_res_with_rownames, "C:/Users/ritas/Desktop/TESE_IMM/GSEA/primary_vs_metastase/GSEA_RES_DEGS.xlsx")

#####################################################################################################33

signficant_res <-as.data.frame(signficant_res[order(signficant_res$padj),])
write.table(signficant_res,"C:/Users/ritas/Desktop/TESE_IMM/cibersort_my_lung_samples/DESEQ2/DESEQ_2/list_significant_NON_PROGRESSIVE_VS_PROGRESSIVE_GENE_NAME.txt",quote = FALSE,sep=",",row.names=FALSE,col.names = FALSE)
View(signficant_res)


high_expressed_genes_2_<-subset(signficant_res,log2FoldChange>0)
high_expressed_genes_2_<-as.data.frame(high_expressed_genes_2_)
View(high_expressed_genes_2_)

low_expressed_genes_2_<-subset(signficant_res,log2FoldChange<0)
low_expressed_genes_2_<-as.data.frame(low_expressed_genes_2_)
View(low_expressed_genes_2_)

##############################################################################
rlog<- rlog(ddsHTSeq__LUNG_filtered)
head(assay(rlog), 3)

#It uses the design formula to calculate the within-group variability (if blind=FALSE) or the across-all-samples variability (if blind=TRUE).
vsd<-vst(ddsHTSeq__LUNG_filtered,blind = FALSE)
head(assay(vsd), 3)
vst <- assay(vst(ddsHTSeq__LUNG_filtered))


library(tidyverse)
library(stats)
#remove batch effect:
mat <- assay(vsd)
mm <- model.matrix(~RESPONSE, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$batch, design=mm)
assay(vsd) <- mat
View(vsd)

######################DO THE PCA##################################
f<-vsd[,vsd@colData@rownames!="Lung9"]#TRY TO REMOVE THIS, TO SEE IF IT HAS RELATIVE CHANGE IN PCA-BIT NO
c<-plotPCA(vsd,intgroup="TYPE")
c


library(ggplot2)
library(PCAtools)

####################DO THE SCREE PLOT############################
p <- pca(vst)
p
screeplot(p,title = "SCREE PLOT 2ºLINE PATIENTS")
pairsplot(p)#mudar cor
pca_result <- prcomp(t(assay(vsd)))


##################################################################################
#PARTE2:VISUALIZATION
#A.MA plot of RNAseq data for entire dataset
#1º Before making the MA-plot, we use the lfcShrink function to shrink the log2 fold changes for the comparison of dex treated vs untreated samples:
library("apeglm")
resultsNames(ddsHTSeq__LUNG_filtered)#[1] "Intercept" "TREATMENT_GROUPS_2_vs_1" "TREATMENT_GROUPS_3_vs_1"
res_2_apeglm<- lfcShrink(ddsHTSeq__LUNG_filtered, coef=2, type="apeglm",lfcThreshold=1)
plotMA(res_2_apeglm, ylim = c(-20, 20))

#B:VOLCANO PLOT:
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
significance_threshold <- 0.05

# Define color vector
keyvals <- ifelse(res$padj > significance_threshold, 'darkgray',
                  ifelse(res$log2FoldChange < 0, 'steelblue', 'tomato'))

keyvals[is.na(keyvals)] <- 'darkgray'
names(keyvals)[keyvals == 'tomato'] <- 'High'
names(keyvals)[keyvals == 'darkgray'] <- 'NS'
names(keyvals)[keyvals == 'steelblue'] <- 'Low'


library(EnhancedVolcano)
EnhancedVolcano(res,
                lab =  ifelse(res$padj < significance_threshold, rownames(res), ""),title = "METASTATIC VS PRIMARY SAMPLES",
                x = 'log2FoldChange',
                y = 'padj',pCutoff = 0.05,colAlpha = 1,colCustom = keyvals,selectLab = rownames(res)[which(names(keyvals) %in% c('High',"NS",'Low'))],ylim=c(0,20),xlim=c(-15,15))

#C.Gene clustering:
library( "genefilter" )
library("RColorBrewer")
library("gplots")
anno <- as.data.frame(colData(vsd)[,"TYPE"])
topVarGenes <- head( order( rowVars( assay(vsd) ), decreasing=TRUE ), 50 )

par(mar=c(10,4,4,2))
par(oma=c(1,15,1,15))
margins = c(8, 8)
lhei=c(10, 10)

# Adjust the log2FoldChange threshold and padj threshold as needed
DEGs <- subset(results(ddsHTSeq__LUNG_filtered), padj < 0.05 & abs(log2FoldChange) > 0.5)
DEGs
View(DEG_names)
# Extract the gene names of DEGs
DEG_names <- rownames(DEGs)

# Replace "vsd" with your actual data
# Assuming "vsd" is a DESeqDataSet object
top_50_DEG_names <- DEG_names[1:30]
heatmap_data <- assay(vsd)[top_50_DEG_names, ]

heatmap.2( heatmap_data, scale="row", 
           trace="none", dendrogram="column",
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           ColSideColors = c( "1"="darkolivegreen3", "2"="maroon")[
             colData(vsd)[,"TYPE"]],margins = c(8, 16),cexRow = 1.2,cexCol = 1.5)


legend("left", inset = c(-0.089, 0.8), legend=c("METASTATIC","PRIMARY"),col = c("darkolivegreen3","maroon"), 
       lty=1, lwd=3, cex=0.6)

legend(y=0.8, x=-0.089, xpd=TRUE,     
       legend=c("METASTASE","PRIMARY"),title = "PRIMARY vs METASTASE",
       col = c("darkolivegreen3","maroon"), 
       lty= 1,             
       lwd = 5,           
       cex=0.5
)
