########## ########## ########## #####
# TXIMPORT AND AGGREGATION OF COUNTS #
########## ########## ########## #####

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("tximport")
#BiocManager::install("DESeq2")
#BiocManager::install("apeglm")
#BiocManager::install("EnhancedVolcano")
#BiocManager::install("edgeR")
#install.packages("RColorBrewer")

library(tximport) #for tximport() function
library(DESeq2) #for DESeqDataSetFromMatrix() function
library(edgeR)  #function "filterByExpr"
library(pheatmap)
library(FactoMineR) #function "PCA"
library(factoextra) #function "fviz_pca_ind"
library(apeglm) #for lfcShrink(x, apeglm...)
library(EnhancedVolcano)  #for volcano plots
library(edgeR) #for filterByExpr() function
library(RColorBrewer)
library(stringr) #for str_split_fixed() function


# Set up the path for the project:
getwd()
# If we are in "scr" folder, move up a folder by running:
setwd("..")
getwd()  ## COMPROBAR QUE ES EL DIRECTORIO QUE QUEREMOS!!!
basedir <- getwd()

# Set up the directory for the results:
resultsDirectory <- file.path(basedir, "PDF/")
tablesDirectory <- file.path(basedir, "Tables/")
#databasesDirectory <- file.path(basedir, "Databases/")
#o
#resultsDirectory <- "/usr/local/Projects/Macroadenomas//PDF"
#tablesDirectory <- "/usr/local/Projects/Macroadenomas/Tables"
          

########## ########## ########## #####
# 1- create a phenodata table and a vector pointing to the quantification files
########## ########## ########## #####

pheno <- read.csv2(file = paste(tablesDirectory, "PhenoData_macroadenomas_v3.csv", sep = "/"), sep = ",")

pheno$Sent_RNA_conc <- as.numeric(pheno$Sent_RNA_conc)
pheno$Valoracion_RNA <- as.numeric(pheno$Valoracion_RNA)
pheno$GW_RNA_Qubit <- as.numeric(pheno$GW_RNA_Qubit)
pheno$GW_RNA_conc <- as.numeric(pheno$GW_RNA_conc)
pheno$GW_RNA_RQN <- as.numeric(pheno$GW_RNA_RQN)
pheno$Ratio_28S.18S <- as.numeric(pheno$Ratio_28S.18S)
pheno$Tamano_del_tumor_mm <- as.numeric(pheno$Tamano_del_tumor_mm)
pheno$Ki67_porcentaje <- as.numeric(pheno$Ki67_porcentaje)

pheno$Fecha_de_nacimiento <- as.Date(pheno$Fecha_de_nacimiento, origin = "1899-12-30") 
pheno$Fecha_de_recurrencia <- as.Date(pheno$Fecha_de_recurrencia, origin = "1899-12-30") 
pheno$Fecha_de_recidiva <- as.Date(pheno$Fecha_de_recidiva, origin = "1899-12-30") 

pheno$Sexo <- as.factor(pheno$Sexo)
pheno$Subclasificacion_histologica <- as.factor(pheno$Subclasificacion_histologica)
pheno$Recurrencia <- as.factor(pheno$Recurrencia)
pheno$Recidiva <- as.factor(pheno$Recidiva)
pheno$Recidiva_O_Recurrencia <- as.factor(pheno$Recidiva_O_Recurrencia)
pheno$Recidiva_Y_Recurrencia <- as.factor(pheno$Recidiva_Y_Recurrencia)
pheno$Grado_de_invasion_escala_Knosp <- as.factor(pheno$Grado_de_invasion_escala_Knosp)
pheno$Invasion_Seno_Cavernoso <- as.factor(pheno$Invasion_Seno_Cavernoso)
pheno$Grado_de_reseccion <- as.factor(pheno$Grado_de_reseccion)
pheno$Fumador <- as.factor(pheno$Fumador)

summary(pheno)

dir <- file.path(basedir, "Genes_results/")
file <- file.path(dir,paste(pheno$Sample,".genes.results", sep =""))       
names(file) <- pheno$Sample


########## ########## ########## #####                
# 2- Use tximport to load gene estimated counts outputs / filter genes
########## ########## ########## #####

# Tximport the gene counts, not "directly" correcting the counts (see $countsFromAbundance == "no")
rsem.genes <- tximport(files = file, type = "rsem", ignoreTxVersion = F, countsFromAbundance = "no", abundanceCol = "TPM")
# We used countsFromAbundance = "no" (default) because we want the original unnormalized gene-level (or transcript-level) counts

# We get a list with estimated counts, effective lengths, abundance. Tximport did not produce warnings, this means our annotation tx2gene matches our salmon files perfectly
head(rsem.genes$counts)

# Tximport to DESeq2
# Create the DESeq2 object with the data, while correcting for isoform length bias  
dds.genes <- DESeqDataSetFromMatrix(countData = round(rsem.genes$counts,0), colData = pheno, design = ~Recidiva_O_Recurrencia)
# We now have our data in DESeq2 format for exploratory and DGE/DTE analyses
# For gene-level counts, the correction is carried "somewhere" inside the object, I believe in the gene length matrix
head(dds.genes@colData)

# Counts
head(rsem.genes$counts)
head(counts(dds.genes))

# Lengths for each sample
head(rsem.genes$length)
#head(assays(dds.genes)[["avgTxLength"]])


########## ########## ########## #####                
# 3- Perform exploratory analysis 
########## ########## ########## #####  

# Prefiltering

# DESeq2 does it internally through the results() function. It filters out:
# genes with low expression: they affect multiple testing
# genes with outlier high expression

# However we will do it here for exploratory analysis, using edgeR's filterByExpr to remove highly and low expressed genes / txs:
cts.genes <- counts(dds.genes)

cts.genes.f <- filterByExpr(y = cts.genes, design = model.matrix(~0+Recidiva_O_Recurrencia, data = dds.genes@colData), group = NULL, lib.size = NULL,
                                      min.count = 10, min.total.count = 200) # This parameters can be changed

dds.genes.f <- dds.genes[as.vector(cts.genes.f),]
cts.genes.f <- cts.genes[as.vector(cts.genes.f),]

# Explore: see that there are some high-count outliers, which we are protected from when using log(counts)
#setwd(resultsDirectory)

pdf(file= paste(resultsDirectory,"1_histograms_counts.pdf",sep="/"), height = 7, width = 5)

layout(matrix(c(1:2),2,1))

hist(cts.genes, breaks = 800000, main = "counts (genes)", xlim= c(0,200), ylim = c(0,1500000),
     cex.lab=0.8, cex.axis=0.8, cex.main=1, cex.sub=0.9)    
hist(cts.genes.f, breaks = 800000, main = "filtered counts (genes)", xlim= c(0,200), ylim = c(0,80000),
     cex.lab=0.8, cex.axis=0.8, cex.main=1, cex.sub=0.9)

dev.off()

par(mfrow=c(1,1))
# MDS plots
pdf(file= paste(resultsDirectory,"1_MDS_plots_counts.pdf",sep="/"), height = 7, width = 7)
# Gene-level counts
# All groups
plotMDS(cts.genes.f, top = nrow(cts.genes.f), gene.selection = "common", col=c("darkgray","darkorange2")[factor(dds.genes@colData$Recidiva_O_Recurrencia)], main = "cts, gene-level")
plotMDS(log(cts.genes.f+1,2), top = nrow(cts.genes.f), gene.selection = "common", col=c("darkgray","darkorange2")[factor(dds.genes@colData$Recidiva_O_Recurrencia)], main = "Log cts, gene-level")

dev.off()
dev.off()

# Heatmaps of correlations between samples
pdf(file= paste(resultsDirectory,"1_heatmap_plots_counts.pdf",sep="/"), height = 2.5, width = 3)
pheatmap::pheatmap(cor(cts.genes),cluster_rows = T, cluster_cols = T, main = "pearson cor, cts, gene-level",fontsize = 3)
dev.off()


########## ########## ########## ########## ##########      
# 4- Perform DTE: intercept model, no subsetting, genes
########## ########## ########## ########## ##########  
# DGE (differential gene expression) by construction, 
# the expressed transcripts were also DTE (differential transcript expression), 
# but they did not count as DTU (differential transcript usage), 
# as the proportions within the gene remained constant. 

# We will split our data for the comparisons.
dds.genes.f.R_vs_No <- dds.genes.f[,c(1:37)]

# Check colData
colData(dds.genes.f.R_vs_No)

# Check design: it is a design with intercept
design(dds.genes.f.R_vs_No)

### Model description.

# Our model is:
# y ~ I0 + b1*x1

# Where:
# I0 is mean of reference group (top factor)
# b1 is difference in mean of group x1 vs reference group

# The order of the factors will determine the reference group used for comparison (the intercept, which is usually the first factor).
# And because we split our data, we may have to RELEVEL the factors to remove unused factors:
colData(dds.genes.f.R_vs_No)$Recidiva_O_Recurrencia; dds.genes.f.R_vs_No$Recidiva_O_Recurrencia # The order of the factor

# Relevel the objects:
dds.genes.f.R_vs_No$Recidiva_O_Recurrencia <- droplevels(dds.genes.f.R_vs_No$Recidiva_O_Recurrencia)
dds.genes.f.R_vs_No$Recidiva_O_Recurrencia <- factor(dds.genes.f.R_vs_No$Recidiva_O_Recurrencia, levels = c('NO','SI'))
#droplevels() function is used to remove unneeded factor levels. This function comes in 
#handy when we need to get rid of factor levels that are no longer in use as a result of 
#subsetting a vector or a data frame.

# Thus, the internally-used design matrix would be
model.matrix(~Recidiva_O_Recurrencia,colData(dds.genes.f.R_vs_No))

### Perform the analyses: DESeq will perform normalization + dispersion + fitting of model
dds.genes.f.R_vs_No <- DESeq(dds.genes.f.R_vs_No, betaPrior = FALSE, test = "Wald")
# betaPrior = FALSE is now the default, and lfcShrink is used afterwards
# test = "Wald" is the default


### Extract results for different comparisons, by using contrast = c(variable, numerator, denominator):
# So if we want: logFC to be, POSITIVE when UP in treatment vs Treat_C, NEGATIVE when DOWN in treatment vs Treat_C,
# it will be: c(group,treatment,Treat_C), because treatment/Treat_C will be > 1 when UP, so logFC > 0, and < 1 when DOWN, so logFC < 0 (logFC(1) is exactly 0)
# independentFiltering of outliers is ON by default. alpha for independent filtering should be the same as the FDR for DGE
res.genes.f.R_vs_No <- results(dds.genes.f.R_vs_No, contrast=c("Recidiva_O_Recurrencia","SI","NO"), independentFiltering = TRUE, alpha = 0.05, lfcThreshold = 0)
# We can do other comparisons using other variables, after setting the model accordingly.

# Summary of results (for each comparison)
summary(res.genes.f.R_vs_No, alpha=0.05)
table(res.genes.f.R_vs_No$padj < 0.05) 

# Export results as data.frame
res.genes.f.R_vs_No.df <- as.data.frame(res.genes.f.R_vs_No)

# As a result we have six columns of information reported for each gene (row). We can use 
# the mcols() function to extract information on what the values stored in each column represent:
mcols(res.genes.f.R_vs_No, use.names=T)

### Also report LFC shrinkage          

# We will use the recommended "apeglm" estimator, but it requires the use of "coef" syntax instead of "contrast"
resultsNames(dds.genes.f.R_vs_No)
# with this code we get the name of the coef= 'Recidiva_O_Recurrencia_SI_vs_NO' that has 
# to be added to see the next function lfcShrink())

res.genes.f.R_vs_No.LFC <- lfcShrink(dds.genes.f.R_vs_No, type="apeglm", res = res.genes.f.R_vs_No, coef = 'Recidiva_O_Recurrencia_SI_vs_NO') 
#apeglm is the recommended estimator, and requires use with "coef" syntax instead of "contrast"
#using 'apeglm' for LFC shrinkage. If used in published research, please cite:
#  Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
#sequence count data: removing the noise and preserving large differences.
#Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895

# Export shrunk results as data.frame
res.genes.f.R_vs_No.LFC.df <- as.data.frame(res.genes.f.R_vs_No.LFC)
res.genes.f.R_vs_No.LFC.df$ID <- rownames(res.genes.f.R_vs_No.LFC.df)
res.genes.f.R_vs_No.LFC.df <- dplyr::arrange(res.genes.f.R_vs_No.LFC.df,padj)

### Export and plot results     

# Complete the DEG table                    
res.genes.f.R_vs_No.df$ID <- rownames(res.genes.f.R_vs_No.df)
res.genes.f.R_vs_No.df <- dplyr::arrange(res.genes.f.R_vs_No.df,padj)


# Check the direction of change
head(res.genes.f.R_vs_No.df)
# En este paso podemos mirar los dataframes res.genes.XX_vs_CTROL.df y seleccionar el elemento
# que cumple "loglogFold is (+)" o "(-)" para comprobar los counts que tiene respecto al CTROL.
# Hemos visto que el elemento [1] tiene lof2FoldChange (+) mientras que el [4] tiene lof2FoldChange (-)
# As? que vamos a comprobarlo:
counts(dds.genes.f.R_vs_No, normalized = T)[res.genes.f.R_vs_No.df$ID[1],] # up (R_vs_No) and logFold is (+)
counts(dds.genes.f.R_vs_No, normalized = T)[res.genes.f.R_vs_No.df$ID[4],] # down (R_vs_No) and logFold is (-)

# Add the variable of direction of change         
res.genes.f.R_vs_No.df$dir <- "up"
res.genes.f.R_vs_No.df$dir[res.genes.f.R_vs_No.df$log2FoldChange < 0] <- "down"

# Add the shrunk logFC variable
table(res.genes.f.R_vs_No.df$ID == res.genes.f.R_vs_No.LFC.df$ID)
res.genes.f.R_vs_No.df$log2FoldChange_shrunk <- res.genes.f.R_vs_No.LFC.df$log2FoldChange

# Write tables:
#setwd(tablesDirectory)
# DGE tables
write.table(res.genes.f.R_vs_No.df, file= paste(tablesDirectory,"1_DESeq_genelevel_R_vs_No.tsv", sep="/"), sep = "\t", quote = F, col.names = T, row.names = F)
write.table(res.genes.f.R_vs_No.LFC.df, file= paste(tablesDirectory,"1_DESeq_genelevel_R_vs_No_shrink.tsv", sep="/"), sep = "\t", quote = F, col.names = T, row.names = F)


##############################################################
# Expression measurements: counts, normalized and not
dds.genes.f2 <- DESeq(dds.genes.f, betaPrior = FALSE, test = "Wald")
write.table(counts(dds.genes.f2, normalized = F), file = paste(tablesDirectory,"1_genes_counts_filtered.tsv", sep="/"), quote = F, col.names = T, row.names = T, sep = "\t" )
write.table(counts(dds.genes.f2, normalized = T), file = paste(tablesDirectory,"1_genes_counts_filtered_DESeq2_normalized.tsv", sep="/"), quote = F, col.names = T, row.names = T, sep = "\t" )

# rlog counts, taking into account experiment design and not
write.table(assay(rlog(dds.genes.f2, blind = TRUE)), file = paste(tablesDirectory,"1_genes_rlog_filtered.tsv", sep="/"), quote = F, col.names = T, row.names = T, sep = "\t" )
write.table(assay(rlog(dds.genes.f2, blind = FALSE)), file = paste(tablesDirectory,"1_genes_rlog_filtered_DESeq2_normalized.tsv", sep="/"), quote = F, col.names = T, row.names = T, sep = "\t" )
##############################################################


# Do some plots
#setwd(resultsDirectory)
pdf(file= paste(resultsDirectory,"2_DESeq_genelevel.pdf", sep = "/"), height = 6, width = 6)

# P-value histograms
hist(res.genes.f.R_vs_No$pvalue, breaks = 100, main = "DGE pvalue histogram, R_vs_No",
     cex.lab=0.8, cex.axis=0.8, cex.main=1, cex.sub=0.9)

# MAplot
DESeq2::plotMA(res.genes.f.R_vs_No, main = "MAplot DGE, R_vs_No",
               cex.lab=0.8, cex.axis=0.8, cex.main=1, cex.sub=0.9)

DESeq2::plotMA(res.genes.f.R_vs_No.LFC, main = "MAplot DGE, R_vs_No (after shrinkage)",
               cex.lab=0.8, cex.axis=0.8, cex.main=1, cex.sub=0.9)

# Look at the shrunk LFCs vs un-shrunk LFCs: for significant DGEs, they hopefully aren't shrunk to 0 ...
table(rownames(res.genes.f.R_vs_No)==rownames(res.genes.f.R_vs_No.LFC))
plot(x = res.genes.f.R_vs_No$log2FoldChange, y = res.genes.f.R_vs_No.LFC$log2FoldChange, pch = 20,
     cex =c(0.75,1)[factor(res.genes.f.R_vs_No$padj < 0.05 & !is.na(res.genes.f.R_vs_No$padj))],
     col =c(rgb(0,0,0,0.25),rgb(1,0,0,0.5))[factor(res.genes.f.R_vs_No$padj < 0.05 & !is.na(res.genes.f.R_vs_No$padj))], main = "", xlab = "", ylab = "")
title(main="DGEs, R_vs_No", xlab = "LFC", ylab = "shrunk LFC", cex.main = 1, cex.lab = 0.8)
abline(a = 0, b = 1, col = "darkblue", lty = 2)

# Plot top 6 significant results
layout(matrix(c(1:6),nrow=3,ncol=2))
# First, we calculate size factors, add normalizationFactors, to avoid set normalized=FALSE
dds.genes.f.plot <- estimateSizeFactors(dds.genes.f)
for(i in 1:6){
  plot(y = counts(dds.genes.f.plot, normalized = T)[res.genes.f.R_vs_No.df$ID[i],c(1:37)], x = as.factor(colnames(dds.genes.f)[c(1:37)]))
  title(main = paste("top R_vs_No_norm DGE,\n #",res.genes.f.R_vs_No.df$ID[i]), ylab = "norm cts", xlab = "")
} 
# I tried changing normalized = FALSE because I got the following error: 
## Error in h(simpleError(msg, call)) : 
## error in evaluating the argument 'y' in selecting a method for function 'plot': 
## first calculate size factors, add normalizationFactors, or set normalized=FALSE

# I have a problem with the examples shown. The huge variability given by one sample in the 
#"No" group of samples makes that the transcript appears as different when it is not.

# Volcano Plot
library(EnhancedVolcano)
EnhancedVolcano(res.genes.f.R_vs_No,
                lab = rownames(res.genes.f.R_vs_No),
                x = 'log2FoldChange',
                y = 'pvalue', 
                title = 'R_vs_No DGE',
                pCutoff = 10e-5,
                FCcutoff = 1.5,
                selectLab = NA,
                labSize = 1.0,
                #legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 3.0,
                axisLabSize = 10,
                titleLabSize = 12) 
                #ylim = c(0, 50),
                #xlim = c(-10, 10))
                
dev.off()  


# Represent PCAs and barplots:
# Write tables:
#setwd(tablesDirectory)
vst_genes <- as.data.frame(assay(vst(dds.genes.f2, blind = TRUE)))
write.table(vst_genes, file=paste(tablesDirectory, "vst_expression_genes_Macroadenomas.tsv", sep = "/"), sep = "\t", quote = F, col.names = T, row.names = T)

vst_genes <- vst_genes[!is.infinite(rowSums(vst_genes)),]
vst_genes <- vst_genes[head(order(-rowVars(as.matrix(vst_genes))), 250),] #han salido 238 DEGS
res.pca.exp <- PCA(t(vst_genes), graph=FALSE)

pdf(file = paste(resultsDirectory, "PCA_exp_R_vs_No.pdf", sep = "/"), height = 5, width = 5)
fviz_pca_ind(res.pca.exp, repel = TRUE, geom.ind = c("point", "text"), fill.ind = factor(c(rep("Tumor", 30), rep("Recid/Recur", 7))), levels = c("Tumor", "Recid/Recur"), col.ind = "black", pointshape = 21, pointsize = 3, addEllipses = F, palette = "jco") + labs(title ="PCA - Exp - 1K - Normalized")
barplot(height = c(nrow(res.genes.f.R_vs_No.df[which(res.genes.f.R_vs_No.df$dir == "up" & res.genes.f.R_vs_No.df$padj < 0.05),]),
                   nrow(res.genes.f.R_vs_No.df[which(res.genes.f.R_vs_No.df$dir == "down" & res.genes.f.R_vs_No.df$padj < 0.05),])))

#barplot(height = c(nrow(res.genes.ClusterA_vs_ClusterB_T.df[which(res.genes.ClusterA_vs_ClusterB_T.df$dir == "up" & res.genes.ClusterA_vs_ClusterB_T.df$padj < 0.05),]),
#                   nrow(res.genes.ClusterA_vs_ClusterB_T.df[which(res.genes.ClusterA_vs_ClusterB_T.df$dir == "down" & res.genes.ClusterA_vs_ClusterB_T.df$padj < 0.05),])))

dev.off()


# CLUSTERING HEATMAP
# Plot the top 4 most differentially methylated CpGs:
vst_genes <- as.data.frame(assay(vst(dds.genes.f2, blind = TRUE)))
vst_genes <- vst_genes[!is.infinite(rowSums(vst_genes)),]
vst_genes <- vst_genes[head(order(-rowVars(as.matrix(vst_genes))), 250),] #han salido 238 DEGS

annotation_R_vs_No <- data.frame(group = as.character(pheno$Recidiva_O_Recurrencia))
annotation_R_vs_No$group[which(annotation_R_vs_No$group == "No")] <- "No"
annotation_R_vs_No$group[which(annotation_R_vs_No$group == "SI")] <- "R"
rownames(annotation_R_vs_No) <- colnames(vst_genes)

DEGs_sign <- as.data.frame(res.genes.f.R_vs_No.df[which(res.genes.f.R_vs_No.df$pvalue < 0.05),])
DEGs_sign[,c("ENSG_name","Gene_ID")] <- str_split_fixed(DEGs_sign$ID, "_", 2) 
DEGs_sign_ordered <- DEGs_sign[order(abs(DEGs_sign$pvalue)),]
#RNAseq6_850k_diff_ordered[,c("ENSG_ID","Extra")] <- str_split_fixed(RNAseq6_850k_diff_ordered$ENSG_name, ".", 2) 
DEGs_sign_ordered <- DEGs_sign_ordered[!duplicated(DEGs_sign_ordered$ENSG_name), ]
DEGs_sign_ordered_adj <- DEGs_sign_ordered[which(DEGs_sign_ordered$padj < 0.05),]
#rownames(DEGs_sign) <- DEGs_sign$ENSG_name                  

pdf(file= file.path(resultsDirectory, "Clustering_heatmap_expr.pdf"), paper="a4r", height=21, width=28, onefile=TRUE)
par(mfrow=c(1,2))
pheatmap(vst_genes, cluster_rows = T, cluster_cols = T, 
         color = colorRampPalette((brewer.pal(n = 7, name = "Reds")))(100), 
         annotation = annotation_R_vs_No, 
         clustering_method = "ward.D2", 
         main = "Heatmap significant DEGs R_vs_No \n (UnAdj.P-Value = 0.05)")
dev.off()


DEGs_expr <- vst_genes[which(rownames(vst_genes) %in% rownames(DEGs_sign_ordered)),]

pdf(file= file.path(resultsDirectory, "Clustering_heatmap_expr_DEGs.pdf"), paper="a4r", height=21, width=28, onefile=TRUE)
par(mfrow=c(1,2))
pheatmap(DEGs_expr, cluster_rows = T, cluster_cols = T, 
         color = colorRampPalette((brewer.pal(n = 7, name = "Reds")))(100), 
         annotation = annotation_R_vs_No, 
         clustering_method = "ward.D2", 
         main = "Heatmap significant DEGs R_vs_No \n (UnAdj.P-Value = 0.05)")
dev.off()


DEGs_expr_adj <- vst_genes[which(rownames(vst_genes) %in% rownames(DEGs_sign_ordered) & DEGs_sign_ordered$padj < 0.05),]

pdf(file= file.path(resultsDirectory, "Clustering_heatmap_expr_DEGs_adj.pdf"), paper="a4r", height=21, width=28, onefile=TRUE)
par(mfrow=c(1,2))
pheatmap(DEGs_expr_adj, cluster_rows = T, cluster_cols = T, 
         color = colorRampPalette((brewer.pal(n = 7, name = "Reds")))(100), 
         annotation = annotation_R_vs_No, 
         clustering_method = "ward.D2", 
         main = "Heatmap significant DEGs R_vs_No \n (UnAdj.P-Value = 0.05)")
dev.off()


write.table(DEGs_sign_ordered, file=paste(tablesDirectory, "DEGs_sign_ordered.tsv", sep = "/"), sep = "\t", quote = F, col.names = T, row.names = F)
write.table(DEGs_sign_ordered_adj, file=paste(tablesDirectory, "DEGs_sign_ordered_adj.tsv", sep = "/"), sep = "\t", quote = F, col.names = T, row.names = F)


write.table(pheno, file=paste(tablesDirectory, "phenoData_Macroadenomas_RNAseq.tsv", sep = "/"), sep = "\t", quote = F, col.names = T, row.names = F)

