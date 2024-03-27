#CLINICAL DATA
library(data.table)
library(readxl)

progressive_clinical<-read_excel("C:/Users/ritas/Desktop/TESE_IMM/ILLUMINA/BD_Lung_MCPP_Illumina_type_biopsie.xlsx")
View(progressive_clinical)

progressive_clinical$`ID LAB` <- gsub("#", "", progressive_clinical$`ID LAB`)

#ID_progressive<-unique(LC_primary_progressive@data$Tumor_Sample_Barcode)
#ID_non_progressive<-unique(LC_primary_non_progresive@data$Tumor_Sample_Barcode)

#TOTAL_ID<-c(ID_progressive,ID_non_progressive)
library(lubridate)
TOTAL_ID<-c("10_LC","36_LC","38_LC","51_LC","52_LC","59_LC","83_LC","33_LC","35_LC","37_LC","40_LC","42_LC","44_LC","47_LC","49_LC","4_LC","54_LC","63_LC","73_LC","79_LC","7_LC","80_LC","82_LC","9LC")
TOTAL_ID<-as.data.frame(TOTAL_ID)
colnames(TOTAL_ID)<-"ID LAB"
TOTAL_ID$`ID LAB` <- gsub("_", "", TOTAL_ID$`ID LAB`)
View(TOTAL_ID)

# Select rows where the values in "column_name" match each value in the list
progressive_clinical<-as.data.frame(progressive_clinical)

merged_df <- merge(progressive_clinical, TOTAL_ID, by = "ID LAB")
View(merged_df)

library(openxlsx)
#change the format of dates:
#A.DATE OF INICIATION OF IMMUNOTHERAPY
number_date_immunotherapy<-as.numeric(merged_df$`Immuno start date`)

datesConverted <- convertToDate(number_date_immunotherapy); datesConverted

merged_df$`Immuno start date` <- datesConverted

#B.DATE OF PROGRESSION
number_date_progression<-as.numeric(merged_df$`Progression date`)

datesConverted <- convertToDate(number_date_progression); datesConverted

merged_df$`Progression date` <- datesConverted

#C.DATE OF LAST-FOLLOW-UP-NÃO EVENTO
merged_df$`Last follow-up` <- gsub("UTC", "", merged_df$`Last follow-up`)


#DO THE DATES DO DAYS
#DO THE EVENT:
#days_to_progression<-round(difftime(merged_df$`Progression date`,merged_df$`Immuno start date`,units = "days"))

merged_df$"Months to Progression"<-round(interval(merged_df$`Immuno start date`, merged_df$`Progression date`) / months(1),2) 
#DO THE NO EVENT:

merged_df$"Months to last follow up"<-round(interval(merged_df$`Immuno start date`,merged_df$`Last follow-up`) / months(1),2)

#CREATE DATAFRAME:
P_vs_NP_data_frame<-cbind(merged_df$`ID LAB`,
                          merged_df$`Best overall response (BOR)`,
                          merged_df$`Months to Progression`,
                          merged_df$`Months to last follow up`)
colnames(P_vs_NP_data_frame)<-c("ID LAB","Response","Months to Progression","Months to last follow up")
P_vs_NP_data_frame<-as.data.frame(P_vs_NP_data_frame)
View(P_vs_NP_data_frame)

#Kapler Meier:
#TIME<-determine the column time with days to progression if exists progression and days to last followup if doesn't exist progression.
Time<-as.numeric(ifelse(P_vs_NP_data_frame$Response=="PD (Progressive disease)",P_vs_NP_data_frame$`Months to Progression`,P_vs_NP_data_frame$`Months to last follow up`))
Time<-as.data.frame(Time)
colnames(Time)<-"Time"
View(Time)

Progression_status<-as.numeric(ifelse(P_vs_NP_data_frame$Response=="PD (Progressive disease)",1,0))
P_vs_NP_data_frame<-cbind(P_vs_NP_data_frame,Time,Progression_status)
View(P_vs_NP_data_frame)

#MUTATIONS IN PROGRESSIVE
#A.MUTATION GNA14
MUTATION_GNA14<-as.numeric(c("1","0","0","1","0","0","0","0","0","1","0","0","1","1","0","1","0","0","0","0","0","0","1","0"))
MUTATION_GNA14<-as.data.frame(MUTATION_GNA14)
colnames(MUTATION_GNA14)<-"MUTATION GNA14"
P_vs_NP_data_frame<-cbind(P_vs_NP_data_frame,MUTATION_GNA14)
View(P_vs_NP_data_frame)

#B.MUTATION_PRELID1
MUTATION_PRELID1<-as.numeric(c("1","0","0","1","0","0","0","1","0","0","0","0","1","1","0","1","1","0","0","0","0","0","1","0"))
MUTATION_PRELID1<-as.data.frame(MUTATION_PRELID1)
colnames(MUTATION_PRELID1)<-"MUTATION PRELID1"
P_vs_NP_data_frame<-cbind(P_vs_NP_data_frame,MUTATION_PRELID1)
View(P_vs_NP_data_frame)

#C.MUTATION GDF7
MUTATION_GDF7<-as.numeric(c("1","0","0","1","0","1","0","0","0","0","0","0","1","0","0","0","0","0","0","0","0","0","1","0"))
MUTATION_GDF7<-as.data.frame(MUTATION_GDF7)
colnames(MUTATION_GDF7)<-"MUTATION GDF7"
P_vs_NP_data_frame<-cbind(P_vs_NP_data_frame,MUTATION_GDF7)
View(P_vs_NP_data_frame)

#D.MUTATION H2AC21
MUTATION_H2AC21<-as.numeric(c("1","0","0","1","0","0","0","0","0","0","0","0","0","1","0","1","0","0","0","0","0","0","1","0"))
MUTATION_H2AC21<-as.data.frame(MUTATION_H2AC21)
colnames(MUTATION_H2AC21)<-"MUTATION H2AC21"
P_vs_NP_data_frame<-cbind(P_vs_NP_data_frame,MUTATION_H2AC21)
View(P_vs_NP_data_frame)

#F.MUTATION TEX36
MUTATION_TEX36<-as.numeric(c("1","0","0","1","0","1","0","0","0","0","0","0","0","1","0","0","0","0","0","0","0","0","1","0"))
MUTATION_TEX36<-as.data.frame(MUTATION_TEX36)
colnames(MUTATION_TEX36)<-"MUTATION TEX36"
P_vs_NP_data_frame<-cbind(P_vs_NP_data_frame,MUTATION_TEX36)
View(P_vs_NP_data_frame)

#save tsv file:
write.table(P_vs_NP_data_frame,"C:/Users/ritas/Desktop/TESE_IMM/ILLUMINA/maf_clinical.tsv",quote = FALSE,row.names = FALSE,col.names = TRUE,sep="\t")

splots <- list()
#PFS-MUTATION GNA14:
library(survival);
library(survminer)

surv_object <- Surv(unlist(P_vs_NP_data_frame$Time),unlist(P_vs_NP_data_frame$Progression_status))
surv_object

surv_fit_groups<- survfit(surv_object ~ unlist(P_vs_NP_data_frame$`MUTATION GNA14`), data=P_vs_NP_data_frame)

splots[[1]]<-ggsurvplot(surv_fit_groups, conf.int = F, title = "GNA14 MUT/WT",ylab="Progression-free survival (PFS)", xlab = "Time (months)",legend.labs=c("WT", "MUT"),palette=c("#A9A9A9","green3"),censor = T, pval=T,pval.coord = c(65.00 * 0.4, 0.3))

#PFS-MUTATION PRELID1:
library(survival);
library(survminer)

surv_object <- Surv(unlist(P_vs_NP_data_frame$Time),unlist(P_vs_NP_data_frame$Progression_status))
surv_object

surv_fit_groups<- survfit(surv_object ~ unlist(P_vs_NP_data_frame$`MUTATION PRELID1`), data=P_vs_NP_data_frame)

splots[[2]]<-ggsurvplot(surv_fit_groups, conf.int = F, title = "PRELID1 MUT/WT",ylab="Progression-free survival (PFS)",xlab = "Time (months)",legend.labs=c("WT", "MUT"),censor = T,palette=c("#A9A9A9","green3"), pval=T,pval.coord = c(65.00 * 0.4, 0.3))

#PFS-MUTATION GDF7:
library(survival);
library(survminer)

surv_object <- Surv(unlist(P_vs_NP_data_frame$Time),unlist(P_vs_NP_data_frame$Progression_status))
surv_object

surv_fit_groups <- survfit(surv_object ~ unlist(P_vs_NP_data_frame$`MUTATION GDF7`), data=P_vs_NP_data_frame)

splots[[3]]<-ggsurvplot(surv_fit_groups, conf.int = F, title = "GDF7 MUT/WT",ylab="Progression-free survival (PFS)",xlab = "Time (months)",legend.labs=c("WT", "MUT"),censor = T, pval=T,palette=c("#A9A9A9","green3"),pval.coord = c(65.00 * 0.4, 0.3))

#PFS-MUTATION H2AC21:
library(survival);
library(survminer)

surv_object <- Surv(unlist(P_vs_NP_data_frame$Time),unlist(P_vs_NP_data_frame$Progression_status))
surv_object

surv_fit_groups <- survfit(surv_object ~ unlist(P_vs_NP_data_frame$`MUTATION H2AC21`), data=P_vs_NP_data_frame)

splots[[4]]<-ggsurvplot(surv_fit_groups, conf.int = F, title = "H2AC21 MUT/WT",ylab="Progression-free survival (PFS)",xlab = "Time (months)",legend.labs=c("WT", "MUT"),censor = T, pval=T,palette=c("#A9A9A9","green3"),pval.coord = c(65.00 * 0.4, 0.3))

#PFS-MUTATION TEX36:
library(survival);
library(survminer)

surv_object <- Surv(unlist(P_vs_NP_data_frame$Time),unlist(P_vs_NP_data_frame$Progression_status))
surv_object

surv_fit_groups <- survfit(surv_object ~ unlist(P_vs_NP_data_frame$`MUTATION TEX36`), data=P_vs_NP_data_frame)

splots[[5]]<-ggsurvplot(surv_fit_groups, conf.int = F, title = "TEX36 MUT/WT",ylab="Progression-free survival (PFS)",xlab = "Time (months)",legend.labs=c("WT", "MUT"),censor = T, pval=T,palette=c("#A9A9A9","green3"),pval.coord = c(65.00 * 0.4, 0.3))


#join all the pfs plots:
x<-arrange_ggsurvplots(splots, print = TRUE,
                    ncol = 3, nrow = 2)
