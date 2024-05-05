rm(list = ls())

library(readr)
BULKDATA <- read_delim("C:/Users/ritas/Desktop/TESE_IMM/cibersort_my_lung_samples/reads_per_gene_star_MCPP-LUNG/BULKDATA.txt", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)
View(BULKDATA)
BULKDATA$Lung27<-NULL#REMOVE THE SAMPLE 27 BECAUSE IS NOT A LUNG SAMPLE
BULK_DATA<- as.data.frame(BULKDATA[!duplicated(BULKDATA$Gene), ])#SELECT THE UNIQUE GENES
BULK_DATA_rownames <- data.frame(BULK_DATA[,-1], row.names=BULK_DATA[,1])
BULK_DATA_rownames[is.na(BULK_DATA_rownames)] <- 0#CONVERT NA VALUES INTO 0
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
metadata<-SAMPLES_GROUPS_TREATMENT_RNA[complete.cases(SAMPLES_GROUPS_TREATMENT_RNA$RESPONSE),]
metadata<-metadata[complete.cases(metadata$TYPE),]
View(metadata)#34 de 44 patients
metadata$RESPONSE[metadata$RESPONSE=="PARTIAL"]<-"NON-PROGRESSIVE"
metadata$RESPONSE[metadata$RESPONSE=="STABLE"]<-"NON-PROGRESSIVE"
#B.Type of biopsie
metadata<-metadata[metadata$TYPE=="PRIMARIO",]

metadata#25 AMOSTRAS PRIMÁRIAS, EM QUE 1 WAS removed after the 9 LUNG, SO ACTUALLY THEY ARE 24 AMOSTRAS

#C.Type of progressive
#1-LESS PROGRESSIVE
#2-HIGH PROGRESSIVE

rownames(metadata)<-metadata$SAMPLES
#put the samples in the same order in the count data as in the col data:
BULK_DATA_rownames_<-BULK_DATA_rownames[,rownames(metadata)]

#1º CONFIRM IF THE ROW NAMES IN BULK_DATA_ROWNAMES (COUNTS_DATA) ARE ALL PRESENT IN SAMPLES_GROUPS_TREATMENT_RNA (COLDATA)
all(colnames(BULK_DATA_rownames_) %in%rownames(metadata))#TRUE

#2º CONFIRM IF THOSE SAMPLE NAMES ARE IN THE SAME ORDER IN BOTH DATAFRAMES:
all(colnames(BULK_DATA_rownames_)==rownames(metadata))

View(BULK_DATA_rownames)
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

BULK_DATA_rownames_df <- BULK_DATA_rownames_[row.names(BULK_DATA_rownames_) %in% coding_genes_list, ]
View(BULK_DATA_rownames_df)#19,187 genes

#put the samples in the same order in the count data as in the col data:
BULK_DATA_rownames_df_<-BULK_DATA_rownames_df[,rownames(metadata)]

#1º CONFIRM IF THE ROW NAMES IN BULK_DATA_ROWNAMES (COUNTS_DATA) ARE ALL PRESENT IN SAMPLES_GROUPS_TREATMENT_RNA (COLDATA)
all(colnames(BULK_DATA_rownames_df_) %in%rownames(metadata))#TRUE

#2º CONFIRM IF THOSE SAMPLE NAMES ARE IN THE SAME ORDER IN BOTH DATAFRAMES:
all(colnames(BULK_DATA_rownames_df_)==rownames(metadata))

######################################################################################
#START OF THE DESEQ2 ANALYSIS:#PRIMARY
#######################################################################
View(BULK_DATA_rownames_df_)
library(DESeq2)
STAR_LUNG_deseq2<- DESeqDataSetFromMatrix(countData=BULK_DATA_rownames_df_,
                                          colData = metadata,
                                          design=~RESPONSE)
STAR_LUNG_deseq2
View(STAR_LUNG_deseq2)

##############################FILTERING BEFORE FITTING MODEL###########3
# Number of genes before filtering:
nrow(STAR_LUNG_deseq2)#19,187 genes with conversion/63368 with ensembl symbol

# Filter
STAR_LUNG_filtered<- STAR_LUNG_deseq2[rowSums(counts(STAR_LUNG_deseq2)) > 5, ]


# Number of genes left after low-count filtering:
nrow(STAR_LUNG_filtered)#18194 genes with conversion/40948 with ensembl symbol

library(DESeq2)
ddsHTSeq__LUNG_filtered <- DESeq(STAR_LUNG_filtered)
ddsHTSeq__LUNG_filtered
View(ddsHTSeq__LUNG_filtered)
names(ddsHTSeq__LUNG_filtered)

normalized_counts <- counts(ddsHTSeq__LUNG_filtered, normalized = TRUE)

dds <- estimateSizeFactors(dds)

#extract the results in the form of a table using the results function:
res<- results(ddsHTSeq__LUNG_filtered)

View(res)
View(res)
head(res)
res

summary(res)
##############################################################################

signficant_res<-subset(res,padj<0.05)#49 significant genes stable vs progressive
summary(signficant_res)

#WRITE THE RES_1_STABLE_PARTIAL_GENE_NAMES WITH ENSEMBLE NAMES.
all_res<-as.data.frame(res)#18,194 genes DEGS and not DEGS
View(all_res)
write.table(all_res, "C:/Users/ritas/Desktop/TESE_IMM/GSEA/PRIMARY_WITH_RESPONSE_ONLY_CODING_GENES/GSEA_RES.tsv", quote = FALSE, sep = "\t", row.names = TRUE)
#order results by padj value (most significant to least)
library(openxlsx)
# Adding row names as a column
signficant_res_dataframe<-as.data.frame(signficant_res)
signficant_res_with_rownames <- cbind(rownames(signficant_res_dataframe), signficant_res_dataframe)
colnames(signficant_res_with_rownames)[1] <- "GeneName"

write.xlsx(signficant_res_with_rownames, "C:/Users/ritas/Desktop/TESE_IMM/GSEA/PRIMARY_WITH_RESPONSE_ONLY_CODING_GENES/GSEA_RES_DEGS.xlsx")



signficant_res <-as.data.frame(signficant_res[order(signficant_res$padj),])
write.table(signficant_res,"C:/Users/ritas/Desktop/TESE_IMM/GSEA/PRIMARY_WITH_RESPONSE_ONLY_CODING_GENES/GSEA_RES_DEGS.txt",quote = FALSE,sep=",")
View(signficant_res)

high_expressed_genes_2_<-subset(signficant_res,log2FoldChange>0)
high_expressed_genes_2_<-as.data.frame(high_expressed_genes_2_)
View(high_expressed_genes_2_)

low_expressed_genes_2_<-subset(signficant_res,log2FoldChange<0)
low_expressed_genes_2_<-as.data.frame(low_expressed_genes_2_)
View(low_expressed_genes_2_)
list_significant_2_STABLE_PROGRESSIVE<-as.data.frame(rownames(signficant_res))
View(list_significant_2_STABLE_PROGRESSIVE)
write.table(list_significant_2_STABLE_PROGRESSIVE,"C:/Users/ritas/Desktop/TESE_IMM/cibersort_my_lung_samples/DESEQ2/DESEQ_2/list_significant_2_STABLE_PROGRESSIVE_GENE_NAME.txt",quote = FALSE,sep=",",row.names=FALSE,col.names = FALSE)

low_expressed_genes_2_STABLE_PROGRESSIVE<-subset(signficant_res,log2FoldChange<0)
View(low_expressed_genes_2_STABLE_PROGRESSIVE)
low_expressed_genes_2_STABLE_PROGRESSIVE<-as.data.frame(rownames(low_expressed_genes_2_STABLE_PROGRESSIVE))
write.table(low_expressed_genes_2_STABLE_PROGRESSIVE,"C:/Users/ritas/Desktop/TESE_IMM/cibersort_my_lung_samples/DESEQ2/DESEQ_2/low_expressed_genes_2_STABLE_PROGRESSIVE_GENE_NAME.txt",quote=FALSE,sep=",",row.names=FALSE,col.names = FALSE)

high_expressed_genes_2_STABLE_PROGRESSIVE<-subset(signficant_res,log2FoldChange>0)
high_expressed_genes_2_STABLE_PROGRESSIVE<-as.data.frame(rownames(high_expressed_genes_2_STABLE_PROGRESSIVE))
View(high_expressed_genes_2_STABLE_PROGRESSIVE)
write.table(high_expressed_genes_2_STABLE_PROGRESSIVE,"C:/Users/ritas/Desktop/TESE_IMM/cibersort_my_lung_samples/DESEQ2/DESEQ_2/high_expressed_genes_2_STABLE_PROGRESSIVE_GENE_NAME.txt",quote=FALSE,sep=",",row.names=FALSE,col.names = FALSE)


##############################################################################
rlog<- rlog(ddsHTSeq__LUNG_filtered)
colnames(head(assay(rlog), 3))
colnames(BULK_DATA_rownames)

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
c<-plotPCA(vsd,intgroup="RESPONSE",ntop=500)
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
EnhancedVolcano(res,title = "P VS NP SAMPLES",
                lab =  ifelse(res$padj < significance_threshold, rownames(res), ""),
                x = 'log2FoldChange',
                y = 'padj',pCutoff = 0.05,colAlpha = 1,colCustom = keyvals,selectLab = rownames(res)[which(names(keyvals) %in% c('High',"NS",'Low'))],ylim=c(0,30),xlim=c(-15,15))

#C.Gene clustering:
library( "genefilter" )
library("RColorBrewer")
library("gplots")
anno <- as.data.frame(colData(vsd)[,"RESPONSE"])
topVarGenes <- head( order( rowVars( assay(vsd) ), decreasing=TRUE ), 30 )

par(mar=c(10,4,4,2))
par(oma=c(1,15,1,15))
margins = c(8, 8)
lhei=c(10, 10)
# Adjust the log2FoldChange threshold and padj threshold as needed
DEGs <- subset(results(ddsHTSeq__LUNG_filtered), padj < 0.05 & abs(log2FoldChange) > 1)
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
             colData(vsd)[,"RESPONSE"]],margins = c(8, 16),cexRow = 1.2,cexCol = 1.5)

legend("left", inset = c(-0.00089, 0.8), legend=c("NON-PROGRESSIVE","PROGRESSIVE"),col = c("darkolivegreen3","maroon"), 
       lty=1, lwd=4, cex=0.6)


#HEATMAP WITH ONLY THE DEGS:
# Assuming you have already installed and loaded the necessary packages
library(gplots)
library(DESeq2)

# Replace "dds" with your actual DESeqDataSet object
# Assuming you have obtained DEGs using results(dds)
# Adjust the log2FoldChange threshold and padj threshold as needed
DEGs <- subset(results(ddsHTSeq__LUNG_filtered), padj < 0.05 & abs(log2FoldChange) > 1)
DEGs
View(DEG_names)
# Extract the gene names of DEGs
DEG_names <- rownames(DEGs)

# Replace "vsd" with your actual data
# Assuming "vsd" is a DESeqDataSet object
heatmap_data <- assay(vsd)[DEG_names, ]

par(mar = c(5.1, 4.1, 4.1, 2.1))# Adjust margins as needed

# Create the heatmap for DEGs
heatmap.2(
  heatmap_data,
  scale = "row",
  trace = "none",
  dendrogram = "none",
  ColSideColors = c( "1"="darkolivegreen3", "2"="maroon")[
    colData(vsd)[,"STAGE"]],
  col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255),
  margins = c(8, 16)
)
par(mar = c(5.1, 4.1, 4.1, 2.1))  # Adjust margins as needed

# Add the legend
legend(
  y = 0.8, x = -0.0089, xpd = TRUE,
  legend = c("STAGE III", "STAGE IV"), title = "STAGE 4 VS 3",
  col = c("darkolivegreen3", "maroon"),
  lty = 1,
  lwd = 5,
  cex = 0.5
)



legend(y=0.8, x=-0.0089, xpd=TRUE,     
       legend=c("STAGE III","STAGE IV"),title = "STAGE 4 VS 3",
       col = c("darkolivegreen3","maroon"), 
       lty= 1,             
       lwd = 5,           
       cex=.5
)









#GENE CLUSTERING WITHOUT HIERARQUICAL CLUSTERING:
# Assuming vsd is your DESeqDataSet object
anno <- as.data.frame(colData(vsd)[,"STAGE"])
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing=TRUE), 7)

# Perform K-means clustering
k <- 2  # Set the number of clusters
kmeans_result <- kmeans(assay(vsd)[topVarGenes, ], centers = k, nstart = 25)
cluster_assignments <- kmeans_result$cluster

# Subset cluster assignments for the topVarGenes
topVarClusterAssignments <- cluster_assignments[topVarGenes]

# Create a color vector for the clusters
cluster_colors <- c("darkolivegreen3", "maroon")[topVarClusterAssignments]

# Plot heatmap with K-means clustering
par(mar=c(10,4,4,2))
par(oma=c(1,15,1,15))
margins = c(8, 8)
lhei=c(10, 10)
library(gplots)
heatmap.2(assay(vsd)[topVarGenes, ], scale="row", 
          trace="none", dendrogram="none",  # Remove column dendrogram
          col = colorRampPalette(rev(brewer.pal(7, "RdBu")))(200),
          margins = c(8, 16))

# Add legend
legend(y=0.8, x=0.089, xpd=TRUE,     
       legend=c("Cluster 1", "Cluster 2"), title = "Clusters",
       col = c("darkolivegreen3", "maroon"), 
       lty= 1,             
       lwd = 5,           
       cex=.5
)


