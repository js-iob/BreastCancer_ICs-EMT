#Author: K.T.Shreya
#Date: 10/12/2021
#Purpose: Estimate co-expressed ion channels and EMT-related genes from rna-seq and microarray datasets of patients with breast cancer (WGCNA)

#Import libraries
library(WGCNA)
library(flashClust)

#options(stringsAsFactors = T)

#Data pre-processing
dat1 = read.csv("ion_GSE42568_NM_042821.csv", sep = '\t', header = TRUE)
dat2 = read.csv("ion_GSE52604_NM_042821.csv", sep = '\t', header = TRUE)
dat3 = read.csv("ion_brca_DEGs_NM_042821.csv", sep = '\t', header = TRUE)

dat1 = t(dat1)
write.table(dat1, file = 'transposed_ion_GSE42568_DEGs_NM_042821.csv', sep = '\t', row.names=TRUE, col.names=FALSE)
dat1 = read.table("transposed_ion_GSE42568_DEGs_NM_042821.csv", sep = '\t', header = TRUE,row.names = 1)
dat1

dat2 = t(dat2)
write.table(dat2, file = 'transposed_ion_GSE52604_DEGs_NM_042821.csv', sep = '\t', row.names=TRUE, col.names=FALSE)
dat2 = read.table("transposed_ion_GSE52604_DEGs_NM_042821.csv", sep = '\t', header = TRUE, row.names = 1)
dat2

dat3 = t(dat3)
write.table(dat3, file = 'transposed_ion_brca_DEGs_NM_042821.csv', sep = '\t', row.names=TRUE, col.names=FALSE)
dat3 = read.table("transposed_ion_brca_DEGs_NM_042821.csv", sep = '\t', header = TRUE, row.names = 1)
dat3

gsg = goodSamplesGenes(dat1, verbose = 3)
gsg
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(dat1)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(dat1)[!gsg$goodSamples], collapse = ", ")))
  dat1 = dat1[gsg$goodSamples, gsg$goodGenes]
}
dim(dat1)

gsg = goodSamplesGenes(dat2, verbose = 3)
gsg
dim(dat2)

gsg = goodSamplesGenes(dat3, verbose = 3)
gsg
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(dat3)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(dat3)[!gsg$goodSamples], collapse = ", ")))
  dat3 = dat3[gsg$goodSamples, gsg$goodGenes]
}
dim(dat3)

#Soft-thresholding parameter
par(mar=c(5.1, 4.1, 4.1, 2.1))
par(mfrow=c(1,2))
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(dat1, powerVector = powers, verbose = 5, networkType = "signed")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence_GSE42568"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,col="red");
abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity_GSE42568"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,col="red")

par(mfrow=c(1,2))
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(dat2, powerVector = powers, verbose = 5, networkType = "signed")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence_GSE52604"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,col="red");
abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity_GSE52604"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,col="red")

par(mfrow=c(1,2))
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(dat3, powerVector = powers, verbose = 5, networkType = "signed")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence_TCGA"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,col="red");
abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity_TCGA"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,col="red")

#Adjacency matrix
softPower1 = 9
adjacency1 = adjacency(dat1, power = softPower1, type = "signed")
dissTOM1   = 1-TOMsimilarity(adjacency1, TOMType="signed")
geneTree1  = flashClust(as.dist(dissTOM1), method="average")

softPower2 = 10
adjacency2 = adjacency(dat2, power = softPower2, type = "signed")
dissTOM2   = 1-TOMsimilarity(adjacency2, TOMType="signed")
geneTree2  = flashClust(as.dist(dissTOM2), method="average")

softPower3 = 6
adjacency3 = adjacency(dat3, power = softPower3, type = "signed")
dissTOM3   = 1-TOMsimilarity(adjacency3, TOMType="signed")
geneTree3  = flashClust(as.dist(dissTOM3), method="average")

#Heirarchical clustering
par(mfrow=c(1,1))
minModuleSize = 20
dynamicMods1 = cutreeDynamic(dendro = geneTree1, distM = dissTOM1,
                             deepSplit = 4, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
dynamicMods1
table(dynamicMods1)
dynamicColors1 = labels2colors(dynamicMods1)
colors1 = table(dynamicColors1)
colors1

plotDendroAndColors(geneTree1, dynamicColors1, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "GSE42568, No.of Samples = 121")

par(mfrow=c(1,1))
minModuleSize = 20
dynamicMods2 = cutreeDynamic(dendro = geneTree2, distM = dissTOM2,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
dynamicMods2
table(dynamicMods2)
dynamicColors2 = labels2colors(dynamicMods2)
colors2 = table(dynamicColors2)
colors2

plotDendroAndColors(geneTree2, dynamicColors2, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "GSE52604, No.of Samples = 55")

par(mfrow=c(1,1))
minModuleSize = 20
dynamicMods3 = cutreeDynamic(dendro = geneTree3, distM = dissTOM3,
                             deepSplit = 1, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
dynamicMods3
table(dynamicMods3)
dynamicColors3 = labels2colors(dynamicMods3)
colors3 = table(dynamicColors3)
colors3

plotDendroAndColors(geneTree3, dynamicColors3, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "TCGA, No.of Samples = 1389")

#Module trait correlation 1
par(mfrow=c(1,3))
traitdata1 = read.csv("GSE42568_traits.csv", sep = '\t',header = TRUE)
patient1 = rownames(dat1)
traitRows1 = match(patient1, traitdata1$Sample_ID)
datTraits1 = traitdata1[traitRows1, -1]
names(datTraits1)
nGenes1 = ncol(dat1)
nSamples1 = nrow(dat1)
MEs0 = moduleEigengenes(dat1, dynamicColors1)$eigengenes
MEs1 = orderMEs(MEs0)
MEs1
moduleTraitCor1 = cor(MEs1, datTraits1, use = "p")
moduleTraitPvalue1 = corPvalueStudent(moduleTraitCor1, nSamples1)
textMatrix1 = paste(signif(moduleTraitCor1, 2), "\n(",signif(moduleTraitPvalue1, 1), ")", sep = "")
dim(textMatrix1) = dim(moduleTraitCor1)
labeledHeatmap(Matrix = moduleTraitCor1,
               xLabels = names(datTraits1),
               yLabels = names(MEs1),
               ySymbols = NULL,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix1,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1), 
               main = paste("Module- Trait relationships - GSE42568"), ylab = "Gene Expression-Based Modules")

#Module trait correlation 2
traitdata2 = read.csv("GSE52604_traits.csv", sep = '\t',header = TRUE)
patient2 = rownames(dat2)
traitRows2 = match(patient2, traitdata2$Sample_ID)
datTraits2 = traitdata2[traitRows2, -1]
names(datTraits2)
nGenes2 = ncol(dat2)
nSamples2 = nrow(dat2)
MEs0 = moduleEigengenes(dat2, dynamicColors2)$eigengenes
MEs2 = orderMEs(MEs0)
MEs2
moduleTraitCor2 = cor(MEs2, datTraits2, use = "p")
moduleTraitPvalue2 = corPvalueStudent(moduleTraitCor2, nSamples2)
textMatrix2 = paste(signif(moduleTraitCor2, 2), "\n(",signif(moduleTraitPvalue2, 1), ")", sep = "")
dim(textMatrix2) = dim(moduleTraitCor2)
labeledHeatmap(Matrix = moduleTraitCor2,
               xLabels = names(datTraits2),
               yLabels = names(MEs2),
               ySymbols = NULL,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix2,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1), 
               main = paste("Module- Trait relationships-GSE52604"), ylab = "Gene Expression-Based Modules")

#Module trait correlation 3
traitdata3 = read.csv("tcga_traits.csv", sep = '\t',header = TRUE)
patient3 = rownames(dat3)
patient3
traitdata3$Sample_id
traitRows3 = match(patient3, traitdata3$Sample_id)
traitRows3
datTraits3 = traitdata3[traitRows3, -1]
names(datTraits3)
nGenes3 = ncol(dat3)
nSamples3 = nrow(dat3)
MEs0 = moduleEigengenes(dat3, dynamicColors3)$eigengenes
MEs3 = orderMEs(MEs0)
MEs3
moduleTraitCor3 = cor(MEs3, datTraits3, use = "p")
moduleTraitPvalue3 = corPvalueStudent(moduleTraitCor3, nSamples3)
textMatrix3 = paste(signif(moduleTraitCor3, 2), "\n(",signif(moduleTraitPvalue3, 1), ")", sep = "")
dim(textMatrix3) = dim(moduleTraitCor3)
labeledHeatmap(Matrix = moduleTraitCor3,
               xLabels = names(datTraits3),
               yLabels = names(MEs3),
               ySymbols = NULL,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix3,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1), 
               main = paste("Module- Trait relationships_TCGA"), ylab = "Gene Expression-Based Modules")

#GS Vs MM plots
par(mfrow = c(1,2))
metastatic1 = as.data.frame(datTraits1$Metastatic)
names(metastatic1) = "metastatic"
modNames1 = substring(names(MEs1), 3)
geneModuleMembership1 = as.data.frame(cor(dat1, MEs1, use = "p"))
MMPvalue1 = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership1), nSamples1))
names(geneModuleMembership1) = paste("MM", modNames1, sep="")
names(MMPvalue1) = paste("p.MM", modNames1, sep="")
geneTraitSignificance1 = as.data.frame(cor(dat1, metastatic1, use = "p"))
GSPvalue1 = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance1), nSamples1))
names(geneTraitSignificance1) = paste("GS.", names(metastatic1), sep="")
names(GSPvalue1) = paste("p.GS.", names(metastatic1), sep="")
module = "yellow"
column = match(module, modNames1)
moduleGenes = dynamicColors1==module
verboseScatterplot(abs(geneModuleMembership1[moduleGenes, column]),
                   abs(geneTraitSignificance1[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Metastatic",
                   main = paste("GSE42568\nModule membership vs. gene significance\n"),
                   cex.main = 0.7, cex.lab = 0.7, cex.axis = 0.7, col = module)

par(mfrow = c(1,2))
metastatic2 = as.data.frame(datTraits2$Metastatic)
names(metastatic2) = "metastatic"
modNames2 = substring(names(MEs2), 3)
geneModuleMembership2 = as.data.frame(cor(dat2, MEs2, use = "p"))
MMPvalue2 = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership2), nSamples2))
names(geneModuleMembership2) = paste("MM", modNames2, sep="")
names(MMPvalue2) = paste("p.MM", modNames2, sep="")
geneTraitSignificance2 = as.data.frame(cor(dat2, metastatic2, use = "p"))
GSPvalue2 = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance2), nSamples2))
names(geneTraitSignificance2) = paste("GS.", names(metastatic2), sep="")
names(GSPvalue2) = paste("p.GS.", names(metastatic2), sep="")
module = "turquoise"
column = match(module, modNames2)
moduleGenes = dynamicColors2==module
verboseScatterplot(abs(geneModuleMembership2[moduleGenes, column]),
                   abs(geneTraitSignificance2[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Metastatic",
                   main = paste("GSE52604\nModule membership vs. gene significance\n"),
                   cex.main = 0.5, cex.lab = 0.5, cex.axis = 0.5, col = module)

par(mfrow = c(1,2))
metastatic3 = as.data.frame(datTraits3$Metastatic)
names(metastatic3) = "metastatic"
modNames3 = substring(names(MEs3), 3)
geneModuleMembership3 = as.data.frame(cor(dat3, MEs3, use = "p"))
MMPvalue3 = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership3), nSamples3))
names(geneModuleMembership3) = paste("MM", modNames3, sep="")
names(MMPvalue3) = paste("p.MM", modNames3, sep="")
geneTraitSignificance3 = as.data.frame(cor(dat3, metastatic3, use = "p"))
GSPvalue3 = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance3), nSamples3))
names(geneTraitSignificance3) = paste("GS.", names(metastatic3), sep="")
names(GSPvalue3) = paste("p.GS.", names(metastatic3), sep="")
module = "blue"
column = match(module, modNames3)
moduleGenes = dynamicColors3==module
verboseScatterplot(abs(geneModuleMembership3[moduleGenes, column]),
                   abs(geneTraitSignificance3[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Metastatic",
                   main = paste("TCGA\nModule membership vs. gene significance\n"),
                   cex.main = 0.5, cex.lab = 0.5, cex.axis = 0.5, col = module)

#Export results to csv files
turquoise1 = names(dat1)[dynamicColors1=="brown"]
yellow2 = names(dat2)[dynamicColors2=="turquoise"]
turquoise3 = names(dat3)[dynamicColors3=="turquoise"]
write.csv(turquoise1, file = "NM_GSE42568_mod.csv")
write.csv(yellow2, file = "NM_GSE52604_mod.csv")
write.csv(turquoise3, file = "NM_tcga_mod.csv")

#Export files for cytoscape input
TOM1 = TOMsimilarityFromExpr(dat1, power = 8)
modules = c("brown", "blue");
genes = names(dat1)
inModule = is.finite(match(moduleColors, modules));
modGenes = genes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modGenes, modGenes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])



