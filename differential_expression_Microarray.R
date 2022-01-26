#Author: K.T.Shreya
#Date: 10/12/2021
#Purpose: Estimate differentially expressed ion channels and EMT-related genes from microarray datasets of patients with breast cancer (limma)

#Importing libraries
#source("https://bioconductor.org/biocLite.R")
#biocLite("limma")
library(limma)
library("ggplot2")
library ("afffy")

#Data import an pre-processing
input_file = "ion_GSE52604_041721.csv"
file = read.table(file = input_file, sep = '\t',header = TRUE)
dim(file)
file_avg = avereps(file, ID = file$probe)
dim(file_avg)
write.table(file_avg, file = 'ion_GSE52604avg_042221.csv', sep = '\t', row.names=FALSE, col.names=TRUE)
file_avg = read.table("ion_GSE52604avg_042221.csv", sep = '\t', header = TRUE, row.names=1, check.names =TRUE)
head(file_avg)
dim(file_avg)

#Design matrix
sample = factor(c(rep("Metastatic",35), rep("Brain",10), rep("Breast",10)))
design.mat = model.matrix(~0+sample)
colnames(design.mat) = levels(sample)
design.mat

#Contrast matrix
contrast.mat = makeContrasts(Diff = ((Metastatic - Brain)-(Metastatic-Breast)), levels = design.mat)
contrast.mat

#Fit Bayes method
fit = lmFit(file_avg, design.mat)
fit2 = contrasts.fit(fit, contrast.mat)
fit3 = eBayes(fit2)
fit3

deg = topTable(fit3,coef = 'Diff', number = nrow(file_avg))

#Upregulated and downregulated gene identification
deg$diffexpressed = "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
deg$diffexpressed[deg$adj.P.Val<0.05 &  deg$logFC > 0.6] = "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
deg$diffexpressed[deg$adj.P.Val<0.05 & deg$logFC < -0.6] = "DOWN"

#Processing for plotting
deg <- cbind(rownames(deg), data.frame(deg, row.names=NULL))
deg
colnames(deg)[1] = "genes"
deg
plot1 = ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val)))+ 
  geom_point(aes(colour=diffexpressed),size=5, alpha=0.7) + scale_colour_manual(values=c("red","grey","green"))+
  labs(title = "Normal Vs Metastatic")+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.5) +
  geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.5)
plot1

#Export results to csv
a = deg[which(deg$diffexpressed=='UP'),]
b = deg[which(deg$diffexpressed=='DOWN'),]
write.table(a, file = 'ion_GSE52604_up_mn_042221.csv', sep = '\t', row.names=FALSE, col.names=TRUE)
write.table(b, file = 'ion_GSE52604_down_mn_042221.csv', sep = '\t', row.names=FALSE, col.names=TRUE)
