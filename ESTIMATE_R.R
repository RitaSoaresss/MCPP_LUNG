library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
help(package="estimate")

estimateScore(metastases_estimate,out.file,platform = "illumina")

metastases_estimate<-fread("C:/Users/ritas/Desktop/TESE_IMM/3_lung_deconvolution/unique_bulk_data_jan_METASTASES.txt")
in.file <- "C:/Users/ritas/Desktop/TESE_IMM/3_lung_deconvolution/unique_bulk_data_jan_PRIMARY.txt"
out.file <- tempfile(pattern="estimate", fileext=".gct")
outputGCT(in.file, out.file)
View(out.file)

out_2.file <- tempfile(pattern="estimate", fileext=".gct")
estimateScore(out.file, out_2.file,platform = "illumina")
LC_estimate<-fread("C:\\Users\\ritas\\AppData\\Local\\Temp\\RtmpqaZowR\\estimate38ac149b6fb9.gct")
View(LC_estimate)

LC_estimate <- LC_estimate[-c(1:2), ]
LC_estimate$V2<-NULL
col_names <- as.character(LC_estimate[1, ])
colnames(LC_estimate)<-col_names
LC_estimate <- LC_estimate[-1, ]
rownames_LC<-as.character(LC_estimate$NAME)
rownames(LC_estimate)<-rownames_LC

col_names <- as.character(LC_estimate[1, ])
LC_estimate_scores<-t(LC_estimate)
colnames(LC_estimate_scores)<-c("StromalScore","ImmuneScore","ESTIMATEScore")
LC_estimate_scores <- LC_estimate_scores[-1, ]
View(LC_estimate_scores)

library(readxl)
SAMPLES_GROUPS_TREATMENT_RNA <- read_excel("C:/Users/ritas/Desktop/TESE_IMM/cibersort_my_lung_samples/SAMPLES_GROUPS_TREATMENT_BIOPSIE_RNA.xlsx")
View(SAMPLES_GROUPS_TREATMENT_RNA)
SAMPLES_GROUPS_TREATMENT_RNA<-as.data.frame(SAMPLES_GROUPS_TREATMENT_RNA)
rownames(SAMPLES_GROUPS_TREATMENT_RNA)<-SAMPLES_GROUPS_TREATMENT_RNA$SAMPLES

#START CREATING THE PLOT FOR ONLY CELL FRATIONS IN METASTASES SAMPLES
metadata<-SAMPLES_GROUPS_TREATMENT_RNA[complete.cases(SAMPLES_GROUPS_TREATMENT_RNA$TYPE),]
#NOTE:IN THE CASE OF PRIMARY SAMPLES
metadata<-metadata[metadata$TYPE=="PRIMARIO",]
metadata<-metadata[complete.cases(metadata$RESPONSE),]
View(metadata)

# Extract the row names from df2
row_names_to_remove <- rownames(metadata)

# Remove rows from df1 where row names match df2
df1_filtered <- LC_estimate_scores[(rownames(LC_estimate_scores) %in% row_names_to_remove),]
View(df1_filtered)

df1_filtered_<-df1_filtered[,rownames(metadata)]

#NOTE METASTASIS VS PRIMARY
df1_filtered<-as.data.frame(df1_filtered)
df1_filtered$TYPE<-c("PRIMARY","METASTASIS","PRIMARY","METASTASIS","PRIMARY","PRIMARY","PRIMARY","METASTASIS","PRIMARY","METASTASIS","PRIMARY","METASTASIS","PRIMARY","PRIMARY","PRIMARY","PRIMARY","PRIMARY","PRIMARY","PRIMARY","PRIMARY","PRIMARY","PRIMARY","PRIMARY","PRIMARY","METASTASIS","PRIMARY","PRIMARY","PRIMARY","PRIMARY","PRIMARY","METASTASIS","PRIMARY","PRIMARY","PRIMARY","PRIMARY","PRIMARY","METASTASIS","PRIMARY","PRIMARY","PRIMARY","PRIMARY","PRIMARY","PRIMARY","PRIMARY")

#A.STROMAL SCORE ESTIMATION

df1_filtered$StromalScore<-as.numeric(df1_filtered$StromalScore)

# Load necessary libraries
library(ggplot2)  # For plotting
library(ggsignif)

#Step 1:CALCULATE THE SIGNIFICANCE LEVEL OF THE DIFFERENCE BETWEEN THE 2 GROUPS OF PATIENTS REGARDING THE TOTAL B CELLS FRACTIONS
test_result_Stromal<- t.test(df1_filtered$StromalScore ~ df1_filtered$TYPE, data = df1_filtered)

group_colors <- c("turquoise4", "orange2")

# Plot ESTIMATE scores per patient using a boxplot or violin plot
STROMAL_plot<-ggplot(df1_filtered, aes(x = df1_filtered$TYPE, y = df1_filtered$StromalScore,fill=df1_filtered$TYPE)) +
  geom_boxplot(show.legend = FALSE) +  # Boxplot
  # geom_violin() +  # Violin plot (you can choose one of these)
  labs(x="",y = "Stromal Score") +
  scale_fill_manual(values = group_colors) +
  geom_signif(comparisons = list(c("METASTASIS", "PRIMARY")), map_signif_level = TRUE, textsize = 4) +
  geom_text(aes(label = paste("p =", signif(test_result_Stromal$p.value, digits = 2))), x = 1.5, y = max(df1_filtered$StromalScore), vjust = -0.5) +
  theme_minimal()

#B.ImmuneScore ESTIMATION

df1_filtered$ImmuneScore<-as.numeric(df1_filtered$ImmuneScore)

# Load necessary libraries
library(ggplot2)  # For plotting
library(ggsignif)

#Step 1:CALCULATE THE SIGNIFICANCE LEVEL OF THE DIFFERENCE BETWEEN THE 2 GROUPS OF PATIENTS REGARDING THE TOTAL B CELLS FRACTIONS
test_result_Immune<- t.test(df1_filtered$ImmuneScore ~ df1_filtered$TYPE, data = df1_filtered)

group_colors <- c("turquoise4", "orange2")

# Plot ESTIMATE scores per patient using a boxplot or violin plot
Immune_plot<-ggplot(df1_filtered, aes(x = df1_filtered$TYPE, y = df1_filtered$ImmuneScore,fill=df1_filtered$TYPE)) +
  geom_boxplot(show.legend = FALSE) +  # Boxplot
  # geom_violin() +  # Violin plot (you can choose one of these)
  labs(x="",y = "Immune Score") +
  scale_fill_manual(values = group_colors) +
  geom_signif(comparisons = list(c("METASTASIS", "PRIMARY")), map_signif_level = TRUE, textsize = 4) +
  geom_text(aes(label = paste("p =", signif(test_result_Immune$p.value, digits = 2))), x = 1.5, y = max(df1_filtered$ImmuneScore), vjust = -0.5) +
  theme_minimal()

#C.ESTIMATEScore ESTIMATION

df1_filtered$ESTIMATEScore<-as.numeric(df1_filtered$ESTIMATEScore)

# Load necessary libraries
library(ggplot2)  # For plotting
library(ggsignif)

#Step 1:CALCULATE THE SIGNIFICANCE LEVEL OF THE DIFFERENCE BETWEEN THE 2 GROUPS OF PATIENTS REGARDING THE TOTAL B CELLS FRACTIONS
test_result_ESTIMATE<- t.test(df1_filtered$ESTIMATEScore ~ df1_filtered$TYPE, data = df1_filtered)

group_colors <- c("turquoise4", "orange2")

# Plot ESTIMATE scores per patient using a boxplot or violin plot
ESTIMATE_plot<-ggplot(df1_filtered, aes(x = df1_filtered$TYPE, y = df1_filtered$ESTIMATEScore,fill=df1_filtered$TYPE)) +
  geom_boxplot(show.legend = FALSE) +  # Boxplot
  # geom_violin() +  # Violin plot (you can choose one of these)
  labs(x="",y = "ESTIMATE Score") +
  scale_fill_manual(values = group_colors) +
  geom_signif(comparisons = list(c("METASTASIS", "PRIMARY")), map_signif_level = TRUE, textsize = 4) +
  geom_text(aes(label = paste("p =", signif(test_result_ESTIMATE$p.value, digits = 2))), x = 1.5, y = max(df1_filtered$ESTIMATEScore), vjust = -0.5) +
  theme_minimal()

library(gridExtra)
combined_plot <- grid.arrange(STROMAL_plot,Immune_plot,ESTIMATE_plot, ncol = 3)

###################################################################################################
######################PRIMARY SAMPLES ANALYSIS ACORDING TO THE RESPONSE ####################################################################
##############################################################################################
in.file <- "C:/Users/ritas/Desktop/TESE_IMM/3_lung_deconvolution/unique_bulk_data_jan_PRIMARY.txt"
out.file <- tempfile(pattern="estimate", fileext=".gct")
outputGCT(in.file, out.file)
View(out.file)

out_2.file <- tempfile(pattern="estimate", fileext=".gct")
estimateScore(out.file, out_2.file,platform = "illumina")
LC_estimate<-fread("C:\\Users\\ritas\\AppData\\Local\\Temp\\RtmpqaZowR\\estimate38ac149b6fb9.gct")
View(LC_estimate)

LC_estimate <- LC_estimate[-c(1:2), ]
LC_estimate$V2<-NULL
col_names <- as.character(LC_estimate[1, ])
colnames(LC_estimate)<-col_names
LC_estimate <- LC_estimate[-1, ]
rownames_LC<-as.character(LC_estimate$NAME)
rownames(LC_estimate)<-rownames_LC

col_names <- as.character(LC_estimate[1, ])
LC_estimate_scores<-t(LC_estimate)
colnames(LC_estimate_scores)<-c("StromalScore","ImmuneScore","ESTIMATEScore")
LC_estimate_scores <- LC_estimate_scores[-1, ]
View(LC_estimate_scores)

library(readxl)
SAMPLES_GROUPS_TREATMENT_RNA <- read_excel("C:/Users/ritas/Desktop/TESE_IMM/cibersort_my_lung_samples/SAMPLES_GROUPS_TREATMENT_BIOPSIE_RNA.xlsx")
View(SAMPLES_GROUPS_TREATMENT_RNA)
SAMPLES_GROUPS_TREATMENT_RNA<-as.data.frame(SAMPLES_GROUPS_TREATMENT_RNA)
rownames(SAMPLES_GROUPS_TREATMENT_RNA)<-SAMPLES_GROUPS_TREATMENT_RNA$SAMPLES

#START CREATING THE PLOT FOR ONLY CELL FRATIONS IN METASTASES SAMPLES
metadata<-SAMPLES_GROUPS_TREATMENT_RNA[complete.cases(SAMPLES_GROUPS_TREATMENT_RNA$TYPE),]
#NOTE:IN THE CASE OF PRIMARY SAMPLES
metadata<-metadata[metadata$TYPE=="PRIMARIO",]
metadata<-metadata[complete.cases(metadata$RESPONSE),]
View(metadata)

# Extract the row names from df2
row_names_to_remove <- rownames(metadata)

# Remove rows from df1 where row names match df2
df1_filtered <- LC_estimate_scores[(rownames(LC_estimate_scores) %in% row_names_to_remove),]
View(df1_filtered)

df1_filtered_<-df1_filtered[,rownames(metadata)]

#NOTE PRIMARY NP AND P
df1_filtered<-as.data.frame(df1_filtered)
df1_filtered$RESPONSE<-c("NP","NP","P","P","NP","NP","P","NP","P","NP","NP","NP","NP","NP","P","P","P","NP","P","P","NP","NP","NP","NP")


#A.STROMAL SCORE ESTIMATION

df1_filtered$StromalScore<-as.numeric(df1_filtered$StromalScore)

# Load necessary libraries
library(ggplot2)  # For plotting
library(ggsignif)

#Step 1:CALCULATE THE SIGNIFICANCE LEVEL OF THE DIFFERENCE BETWEEN THE 2 GROUPS OF PATIENTS REGARDING THE TOTAL B CELLS FRACTIONS
test_result_Stromal<- t.test(df1_filtered$StromalScore ~ df1_filtered$RESPONSE, data = df1_filtered)

group_colors <- c("turquoise4", "orange2")

# Plot ESTIMATE scores per patient using a boxplot or violin plot
STROMAL_plot_RESPONSE<-ggplot(df1_filtered, aes(x = df1_filtered$RESPONSE, y = df1_filtered$StromalScore,fill=df1_filtered$RESPONSE)) +
  geom_boxplot(show.legend = FALSE) +  # Boxplot
  # geom_violin() +  # Violin plot (you can choose one of these)
  labs(x="",y = "Stromal Score") +
  scale_fill_manual(values = group_colors) +
  geom_signif(comparisons = list(c("NP", "P")), map_signif_level = TRUE, textsize = 4) +
  geom_text(aes(label = paste("p =", signif(test_result_Stromal$p.value, digits = 2))), x = 1.5, y = max(df1_filtered$StromalScore), vjust = -0.5) +
  theme_minimal()

#B.ImmuneScore ESTIMATION

df1_filtered$ImmuneScore<-as.numeric(df1_filtered$ImmuneScore)

# Load necessary libraries
library(ggplot2)  # For plotting
library(ggsignif)

#Step 1:CALCULATE THE SIGNIFICANCE LEVEL OF THE DIFFERENCE BETWEEN THE 2 GROUPS OF PATIENTS REGARDING THE TOTAL B CELLS FRACTIONS
test_result_Immune<- t.test(df1_filtered$ImmuneScore ~ df1_filtered$RESPONSE, data = df1_filtered)

group_colors <- c("turquoise4", "orange2")

# Plot ESTIMATE scores per patient using a boxplot or violin plot
Immune_plot_RESPONSE<-ggplot(df1_filtered, aes(x = df1_filtered$RESPONSE, y = df1_filtered$ImmuneScore,fill=df1_filtered$RESPONSE)) +
  geom_boxplot(show.legend = FALSE) +  # Boxplot
  # geom_violin() +  # Violin plot (you can choose one of these)
  labs(x="",y = "Immune Score") +
  scale_fill_manual(values = group_colors) +
  geom_signif(comparisons = list(c("NP", "P")), map_signif_level = TRUE, textsize = 4) +
  geom_text(aes(label = paste("p =", signif(test_result_Immune$p.value, digits = 2))), x = 1.5, y = max(df1_filtered$ImmuneScore), vjust = -0.5) +
  theme_minimal()

#C.ESTIMATEScore ESTIMATION

df1_filtered$ESTIMATEScore<-as.numeric(df1_filtered$ESTIMATEScore)

# Load necessary libraries
library(ggplot2)  # For plotting
library(ggsignif)

#Step 1:CALCULATE THE SIGNIFICANCE LEVEL OF THE DIFFERENCE BETWEEN THE 2 GROUPS OF PATIENTS REGARDING THE TOTAL B CELLS FRACTIONS
test_result_ESTIMATE<- t.test(df1_filtered$ESTIMATEScore ~ df1_filtered$RESPONSE, data = df1_filtered)

group_colors <- c("turquoise4", "orange2")

# Plot ESTIMATE scores per patient using a boxplot or violin plot
ESTIMATE_plot_RESPONSE<-ggplot(df1_filtered, aes(x = df1_filtered$RESPONSE, y = df1_filtered$ESTIMATEScore,fill=df1_filtered$RESPONSE)) +
  geom_boxplot(show.legend = FALSE) +  # Boxplot
  # geom_violin() +  # Violin plot (you can choose one of these)
  labs(x="",y = "ESTIMATE Score") +
  scale_fill_manual(values = group_colors) +
  geom_signif(comparisons = list(c("NP", "P")), map_signif_level = TRUE, textsize = 4) +
  geom_text(aes(label = paste("p =", signif(test_result_ESTIMATE$p.value, digits = 2))), x = 1.5, y = max(df1_filtered$ESTIMATEScore), vjust = -0.5) +
  theme_minimal()

library(gridExtra)
combined_plot <- grid.arrange(STROMAL_plot_RESPONSE,Immune_plot_RESPONSE,ESTIMATE_plot_RESPONSE, ncol = 3)

