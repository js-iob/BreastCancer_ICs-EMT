#Author: K.T.Shreya
#Date: 01/20/2022
#Purpose: 5 year survival analysis of patients with breast cancer using tcga breast cancer ion channel mRNA expression profiles

#setwd('C:\\Users\\91905\\.......\\survival_012022')
#Importing libraries
library(RTCGA.clinical)
library(survival)
library(dplyr)
library(survminer)
library(ggplot2)

#Data import and processing
dim(BRCA.clinical)

rna <- read.table('BRCA_ions.txt', header=T,sep='\t', fill = TRUE)
dim(rna)
rna = distinct(rna)
dim(rna)
colnames(rna)[1] = 'bcr_patient_barcode' #Rename column name with patient ids as the one similar in clinical data file
head(rna)

clin = survivalTCGA(BRCA.clinical)
head(clin)
table(clin$patient.vital_status)

#Combining gene expression and clinical data
clin = rna %>%
  as_tibble() %>%
  select(bcr_patient_barcode, TRPC6, GJA1, GRIN1, CLCN6, SCN3A, SCN2A, TRPC1, CACNA1B, CLIC2,
         GJA4, VDAC2, VDAC3, VDAC1, CLCN3, GJA3, KCNJ11, KCNH2, GRIN2D, KCNJ10, CACNA2D2,CFTR,
         CLIC1, AQP4, AQP5, CACNA1C, KCNN4, ANO1, ITPR3, GJB2, KCNT1, SCNN1G, KCNJ6, HTR3C,
         CLCN6, GLRB, SCN3A,ANO3, ANO5, ANO6, LRRC8C, TRPV3, GABRG3, KCNT1, GJA5, GRIK5, CLCNKB,
         KCNN1) %>%
  mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) %>%
  inner_join(clin, by="bcr_patient_barcode")
names(clin)

dim(clin)
table(clin$patient.vital_status)

#5yr criteria:
all_clin = within(clin, patient.vital_status[patient.vital_status == 1 & times > 1825] <- 0)
all_clin = within(clin, times[times > 1825] <- 1825)
table(all_clin$patient.vital_status)
range (all_clin$times)

#median(all_clin$CLIC1)
#hist(clin$CLIC2) #To check for normal distribution

#Grouping into high and low expression based on median value of expession
gene = cut(all_clin$KCNN1, breaks = c(0, median(all_clin$KCNN1), Inf), labels = c("low", "high"))

#Survival anlysis with specific ion channel
sfit = survfit(Surv(times, patient.vital_status)~gene, data = all_clin)
sfit

#Cox regression
fit = coxph(Surv(times, patient.vital_status)~KCNN1, data = all_clin)
fit

#Plotting
ggsurv = ggsurvplot(sfit, legend.labs=c("Low","High"), legend.title=
             'Expression', title = 'KCNN1', pval = TRUE, pval.method = TRUE, xlab = "Time(in days)" )

ggsurv$plot + ggplot2::annotate("text",x = Inf, y = Inf, vjust = 1, hjust=1, 
                                label = "HR = 0.96 \n p(HR) = 0.61")
