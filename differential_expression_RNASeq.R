#Author: K.T.Shreya
#Date: 10/12/2021
#Purpose: Estimate differentially expressed ion channels and EMT-related genes from rna-seq datasets of patients with breast cancer (DESeq2)

#Import libraries
library("DESeq2")
library("ggplot2")

#Data import and preprocessing
file = read.table("ion_genes_breast_cancer_2.csv", sep = '\t',header = TRUE)
file = file[,4:ncol(file)]
transpose = t(file)
write.table(transpose, file = 'ion_genes_breast_cancer_transposed_042221.csv', sep = '\t', row.names=TRUE, col.names=FALSE)
file2 = read.table("ion_genes_breast_cancer_transposed_042221.csv", sep = '\t', header = TRUE, row.names=1, check.names =TRUE)
dim(file2)

#Data transformation to get raw reads required for DESeq2
file2 = round(2**file2)-1
dim(file2)

#DESeq2
countdata = as.matrix(file2)
condition = factor(c(rep("normal",292), rep("tumor",1090), rep("metastatic",7)))
coldata = data.frame(row.names = colnames(countdata), condition)
ddsFull = DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design =~condition)
dds = DESeq(ddsFull)
res <- results( dds )
summary(res)
res_ordered = res[order(res$padj),]

res1 = results(dds, contrast = c("condition", "tumor", "normal"))
res1_ordered = res1[order(res1$padj),]
res2 = results(dds, contrast = c("condition", "metastatic", "normal"))
res2_ordered = res2[order(res2$padj),]
summary(res1)
summary(res2)

#Upregulated and downregulated identification
res_d = as.data.frame(res_ordered)
res1_d = as.data.frame(res1_ordered)
res2_d = as.data.frame(res2_ordered)

#In tumor Vs metastatic
res_d$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
res_d$diffexpressed[res_d$log2FoldChange > 0.6 & res_d$padj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
res_d$diffexpressed[res_d$log2FoldChange < -0.6 & res_d$pvalue < 0.05] <- "DOWN"

#Processing for plotting
res_d <- cbind(rownames(res_d), data.frame(res_d, row.names=NULL))
colnames(res_d)[1] <- "genes"
plot1 = ggplot(res_d, aes(x = log2FoldChange, y = -log10(padj)))+ 
  geom_point(aes(colour=diffexpressed),size=2, alpha=1) + scale_colour_manual(values=c("red","blue","green"))+
  labs(title = "Tumor_Metastatic")+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.8)
plot1

#Export results to csv files 
a = res_d[which(res_d$diffexpressed=='UP'),]
b = res_d[which(res_d$diffexpressed=='DOWN'),]
write.table(a, file = 'tm_up.csv', sep = '\t', row.names=FALSE, col.names=TRUE)
write.table(b, file = 'tm_down.csv', sep = '\t', row.names=FALSE, col.names=TRUE)

#In normal Vs tumor
res1_d$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
res1_d$diffexpressed[res1_d$log2FoldChange > 0.6 & res1_d$padj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
res1_d$diffexpressed[res1_d$log2FoldChange < -0.6 & res1_d$pvalue < 0.05] <- "DOWN"

#Processing for plotting
res1_d <- cbind(rownames(res1_d), data.frame(res1_d, row.names=NULL))
#res1_d = res1_d[-c(1)]
res1_d
colnames(res1_d)[1] <- "genes"
#res1_d$delabel <- NA
#res1_d$delabel[res1_d$diffexpressed != "NO"] <- res1_d$genes[res1_d$diffexpressed != "NO"]
res1_d
plot2 = ggplot(res1_d, aes(x = log2FoldChange, y = -log10(padj)))+ 
  geom_point(aes(colour=diffexpressed),size=2, alpha=1) + scale_colour_manual(values=c("red","blue","green"))+
  labs(title = "Normal_Tumor")+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.8)
plot2

#Export results to csv files 
c = res1_d[which(res1_d$diffexpressed=='UP'),]
d = res1_d[which(res1_d$diffexpressed=='DOWN'),]
write.table(c, file = 'tn_up.csv', sep = '\t', row.names=FALSE, col.names=TRUE)
write.table(d, file = 'tn_down.csv', sep = '\t', row.names=FALSE, col.names=TRUE)

#In normal Vs metastatic
res2_d$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
res2_d$diffexpressed[res2_d$log2FoldChange > 0.6 & res2_d$padj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
res2_d$diffexpressed[res2_d$log2FoldChange < -0.6 & res2_d$pvalue < 0.05] <- "DOWN"

#Processing for plotting
res2_d <- cbind(rownames(res2_d), data.frame(res2_d, row.names=NULL))
res2_d
colnames(res2_d)[1] <- "genes"
plot3 = ggplot(res2_d, aes(x = log2FoldChange, y = -log10(padj)))+ 
  geom_point(aes(colour=diffexpressed),size=2, alpha=1) + scale_colour_manual(values=c("red","blue","green"))+
  labs(title = "Normal_Metastatic")+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.8)
plot3

#Export results to csv files 
e = res2_d[which(res2_d$diffexpressed=='UP'),]
f = res2_d[which(res2_d$diffexpressed=='DOWN'),]
write.table(e, file = 'nm_up.csv', sep = '\t', row.names=FALSE, col.names=TRUE)
write.table(f, file = 'nm_down.csv', sep = '\t', row.names=FALSE, col.names=TRUE)

