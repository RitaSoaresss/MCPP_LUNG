#CLINICAL DATA
library(data.table)
library(readxl)

primary_clinical<-read_excel("C:/Users/ritas/Desktop/TESE_IMM/ILLUMINA/BD_Lung_MCPP_Illumina_type_biopsie.xlsx")
View(primary_clinical)

primary_clinical<-BD_Lung_MCPP_Illumina_type_biopsie

progressive_clinical$`ID LAB` <- gsub("#", "", progressive_clinical$`ID LAB`)

#ID_progressive<-unique(LC_primary_progressive@data$Tumor_Sample_Barcode)
#ID_non_progressive<-unique(LC_primary_non_progresive@data$Tumor_Sample_Barcode)

#TOTAL_ID<-c(ID_progressive,ID_non_progressive)
library(lubridate)
TOTAL_ID<-c("52_LC","37_LC","51_LC","83_LC","9_LC","38_LC","36_LC","57_LC","49_LC","82_LC","47_LC","32_LC","35_LC","64_LC","42_LC","54_LC","81_LC","59_LC","4_LC","79_LC","44_LC","10_LC","45_LC","78LC","63_LC","33_LC","7_LC","73_LC","68_LC","40_LC","80_LC")
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
#A.DATE OF DIAGNOSIS
datesConverted <- convertToDate(merged_df$`Date of diagnosis`); datesConverted

merged_df$`Date of diagnosis` <- datesConverted

#B.DATE OF DEATH
merged_df$`Death date`<-gsub("UTC", "", merged_df$`Death date`)

#C.DATE OF LAST-FOLLOW-UP-NÃO EVENTO
merged_df$`Last follow-up` <- gsub("UTC", "", merged_df$`Last follow-up`)


#DO THE DATES DO DAYS
#DO THE EVENT:
#days_to_progression<-round(difftime(merged_df$`Progression date`,merged_df$`Immuno start date`,units = "days"))

merged_df$"Months to Death"<-round(interval(merged_df$`Date of diagnosis`, merged_df$`Death date`) / months(1),2) 

#DO THE NO EVENT:
merged_df$"Months to last follow up"<-round(interval(merged_df$`Date of diagnosis`,merged_df$`Last follow-up`) / months(1),2)

#CREATE DATAFRAME:
primary_data_frame<-cbind(merged_df$`ID LAB`,
                          merged_df$Death,
                          merged_df$`Months to Death`,
                          merged_df$`Months to last follow up`)
colnames(primary_data_frame)<-c("ID LAB","Death","Months to Death","Months to last follow up")
primary_data_frame<-as.data.frame(primary_data_frame)
View(primary_data_frame)

primary_data_frame$Death[2]<-"Yes"
primary_data_frame$Death[11]<-"No"

#Kapler Meier:
#TIME<-determine the column time with days to progression if exists progression and days to last followup if doesn't exist progression.
Time<-as.numeric(ifelse(primary_data_frame$Death=="Yes",primary_data_frame$`Months to Death`,primary_data_frame$`Months to last follow up`))
Time<-as.data.frame(Time)
colnames(Time)<-"Time"
View(Time)

Death_status<-as.numeric(ifelse(merged_df$Death=="Yes",1,0))
Death_status<-as.data.frame(Death_status)
primary_data_frame<-cbind(primary_data_frame,Time,Death_status)
View(primary_data_frame)

primary_data_frame$Death_status[2]<-"1"
primary_data_frame$Death_status[11]<-"0"

#MUTATIONS IN PROGRESSIVE
#A.MUTATION AKT1
MUTATION_AKT1<-as.numeric(c("1","1","0","1","1","0","0","0","0","0","0","1","0","0","1","1","0","0","0","0","0","0","0","0","0","0","1","0","0","0","1"))
MUTATION_AKT1<-as.data.frame(MUTATION_AKT1)
colnames(MUTATION_AKT1)<-"MUTATION AKT1"
primary_data_frame<-cbind(primary_data_frame,MUTATION_AKT1)
View(primary_data_frame)

#save tsv file:
write.table(P_vs_NP_data_frame,"C:/Users/ritas/Desktop/TESE_IMM/ILLUMINA/maf_clinical.tsv",quote = FALSE,row.names = FALSE,col.names = TRUE,sep="\t")

#PFS-MUTATION AKT1:
library(survival);
library(survminer)

surv_object <- Surv(unlist(as.numeric(primary_data_frame$Time)),unlist(as.numeric(primary_data_frame$Death_status)))
surv_object

surv_fit_groups<- survfit(surv_object ~ unlist(as.numeric(primary_data_frame$`MUTATION AKT1`)), data=primary_data_frame)

mutation_akt1_os<-ggsurvplot(surv_fit_groups, conf.int = F, title = "AKT1 MUT/WT",ylab="Overall Survival", xlab = "Time (months)",legend.labs=c("WT", "MUT-Del"),censor = T, pval=T,pval.coord = c(95.00 * 0.2, 0.2),palette = c("#757575","darkorange2"))
