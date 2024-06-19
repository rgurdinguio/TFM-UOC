#############################################################################
### PROJECT: Hypophysis Macroadenomas
###
###
### Script for the analysis of Methylation array data (850K)
###          
###
### Author: Rocio G. Urdinguio 
###         
#############################################################################

# Load the required packages:
library("quadprog")
library("minfi")
library("IlluminaHumanMethylationEPICmanifest")
library("RColorBrewer")
library("missMethyl")
library("matrixStats")
library("stringr")
library("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
library("genefilter")
library("ChAMP")
library("gplots")
library("cluster")
library("gridExtra")
library("reshape2")
library("beanplot")
library("lumi")
library("factoextra")

options("width"=110) #  using R in terminal mode on Linux, making easier the inspection of data sets
###################################################################
#####                       READING DATA                      #####
###################################################################

# Set up the path for the project and the IDATs location:
# basedir <- "/media/rauldiul/Expansion/1_ROCIOyJAVI/RG_Macroadenomas"
getwd()
# If we are in "scr" folder, move up a folder by running:
setwd("..")
getwd()
basedir <- getwd()

# Set up the directory for the database filtering tables (multimapping, crossreactive, SNPs...):
dataDirectory <- file.path(basedir, "Databases/Methylation_filters/")

# Set up the directory for the results:
resultsDirectory <- file.path(basedir, "PDF/")
tablesDirectory <- file.path(basedir, "Tables/")

# Read the sample sheet for the experiment:
setwd(basedir)
targets <- read.metharray.sheet("Samplesheets", pattern = ".csv$")

targets_28 <- targets[!(targets$Condition == "11_No" | targets$Condition == "29_No"),]

# Add the information of the slide array for further filtering purposes:
targets_28$SlideArray <- paste(targets_28$Slide, targets_28$Array, sep = "_")
targets_28$Basename <- paste(basedir, "Raw_data_macro", targets_28$SlideArray, sep = "/")
#targets$color <- c("darkslategray", "#e8baba", "#e36d6d", "#7a0804", "#bfd9a7", "#97d162", "#428505", "#b6c3e3", "#6789db", "#063094", "#d9a0d6", "#c458bf", "#6e0b69")

# Read the raw data from the IDAT files:
RGSet_28 <- read.metharray.exp(targets = targets_28, force = T)

# Generate the phenoData_28 for the desired dataset:
phenoData_28 <- pData(RGSet_28)

# Obtain probe design information from the dataset:
manifest_28 <- getManifest(RGSet_28)

# Get annotation from RGSet_28 (EPIC):
annotationEPIC_28 <- getAnnotation(RGSet_28)

# Set the color palette for visual representations:
groups_28 <- factor(phenoData_28$Sample_Group #, levels = c("Control", "S2", "Stip", "N2", "Ntip")
                 )
pal_28 <- factor(phenoData_28$color #, levels = c("darkslategray", "#7a0804", "#428505", "#063094", "#6e0b69")
              )

###################################################################
#####                    QUALITY CONTROL 1                    #####
###################################################################

# Calculate the detection p-values:
detP_28 <- detectionP(RGSet_28)

# Examine mean detection p-values across all samples to identify any failed samples:
# "las=2" is for making label text perpendicular to axis
pdf(file= paste(resultsDirectory, "01-Qual_pvalue_plot_28.pdf", sep = "/"), paper="a4r", height=21, width=28, onefile=TRUE)
par(mfrow=c(1,1))
barplot(colMeans(detP_28), col=phenoData_28$color, las=2, cex.axis=0.8, col.axis="grey40", cex.names=0.5, ylab="Mean detection p-values", names.arg=phenoData_28$Sample_Name, border=phenoData_28$color, space=1.5)
abline(h=0.01, col="green4")
legend("topright", legend=levels(groups_28), fill=levels(pal_28), bg="white")
barplot(colMeans(detP_28), col=phenoData_28$color, las=2, cex.axis=0.8, col.axis="grey40", cex.names=0.5, ylim=c(0,0.002), ylab="Mean detection p-values", names.arg=phenoData_28$Sample_Name, border=phenoData_28$color, space=1.5)
legend("topright", legend=levels(groups_28), fill=levels(pal_28), bg="white")
dev.off()

# Generate the QC report for the dataset:
qcReport(RGSet_28, sampNames=phenoData_28$Sample_Name, sampGroups = phenoData_28$Sample_Group, pdf= paste(resultsDirectory, "qcReport_28.pdf", sep = "/"))

###################################################################
#####                    SAMPLE FILTERING                     #####
###################################################################

# Remove poor quality samples (i.e detP_28 > 0.01 or detP_28 > 0.05) and IVD:
keep <- colMeans(detP_28) < 0.05
RGSet_28_filtered <- RGSet_28[, keep]

# Remove poor quality samples from RGSet_28, targets data, detP_28 and phenoData_28:
targets_28_filtered <- targets_28[keep,]
detP_28_filtered <- detP_28[,keep]
phenoData_28_filtered <- phenoData_28[keep,]

###################################################################
#####         DATA NORMALIZATION AND QUALITY CONTROL 2        #####
###################################################################

# Generate a MethylSet object with the methylated and the unmethylated signals (without normalization):
MSetraw_28 <- preprocessRaw(RGSet_28_filtered)

# Extract and plot quality control information from the MethylSet objects:
qcraw_28 <- getQC(MSetraw_28)

pdf(file= paste(resultsDirectory, "02-plotQC_28.pdf", sep = "/"), paper="a4r", height=21, width=28, onefile=TRUE)
par(mfrow=c(1,1))
badSampleCutoff = 9
plotQC(qcraw_28, badSampleCutoff)
dev.off()

# Generate a MethylSet object with the methylated and the unmethylated signals (normalized):  
MSetsq_28 <- preprocessNoob(RGSet_28_filtered, offset = 15, dyeCorr = TRUE, dyeMethod =  "single")
#MSetsq <- preprocessFunnorm(RGSet_28_filtered, bgCorr = T)
#MSetsq <- preprocessQuantile(RGSet_28_filtered, fixOutliers = T, quantileNormalize = T, sex = c("M", "M", "M", "M", "M", "M"))

# Visualize the data after normalization:
pdf(file= paste(resultsDirectory, "03-norm_step_28.pdf", sep = "/"), paper="a4r", height=21, width=28, onefile=TRUE)
par(mfrow=c(1,2))
densityBeanPlot(getBeta(MSetsq_28), sampGroups=phenoData_28_filtered$Sample_Group, sampNames = phenoData_28_filtered$Sample_Name, main = "Beta Noob norm data", pal = phenoData_28_filtered$color) 

# Visualize the density M Plots getM(MSetsq_28)
par(mfrow=c(1,2))
densityPlot(getBeta(MSetsq_28), sampGroups=phenoData_28_filtered$Sample_Name, main="Noob norm - Beta values",  legend=FALSE, pal=phenoData_28_filtered$color)
legend("topleft", legend=levels(groups_28), fill=levels(pal_28), bg="white")
densityPlot(getM(MSetsq_28), sampGroups=phenoData_28_filtered$Sample_Name, main="Noob norm - M values", legend=FALSE, pal=phenoData_28_filtered$color)
legend("topleft", legend=levels(groups_28), fill=levels(pal_28), bg="white")

dev.off()


###################################################################
#####                    PROBE FILTERING                      #####
###################################################################

# Remove any probes that have failed in more than 1% samples to avoid excluding CpGs because of our large sample
keep2 <- rowSums(detP_28_filtered < 0.01) >= 0.99*ncol (MSetsq_28)
MSetsq_28Flt <- MSetsq_28[keep2,]

# If the data includes males and females, remove probes on sex chromosomes:
keep3 <- !(featureNames(MSetsq_28Flt) %in% annotationEPIC_28$Name[annotationEPIC_28$chr %in% c("chrX", "chrY")])
MSetsq_28Flt <- MSetsq_28Flt[keep3,]

# Remove probes with SNPs at CpG sites:  Obtained from https://github.com/sirselim/illumina450k_filtering
xProbesinSNPs <- read.csv(file=paste(dataDirectory, " Probes_overlapping_genetic_variants_EPIC.csv", sep="/"), stringsAsFactors=FALSE)  
keep5 <- !(featureNames(MSetsq_28Flt) %in% xProbesinSNPs$PROBE)   # Watch out!!, download and include annotation for EPIC array
MSetsq_28Flt <- MSetsq_28Flt[keep5,]

# Exclude cross reactive probes:  Obtained from https://github.com/sirselim/illumina450k_filtering
xReactiveProbes <- read.csv(file=paste(dataDirectory, " Cross-reactive_probes_on_EPIC_array.csv", sep="/"), stringsAsFactors=FALSE)  
keep6 <- !(featureNames(MSetsq_28Flt) %in% xReactiveProbes$X)   
MSetsq_28Flt <- MSetsq_28Flt[keep6,]

# Exclude multi-mapping probes:  450K, Obtained from https://github.com/sirselim/illumina450k_filtering
multimapProbes <- read.csv(file=paste(dataDirectory, "HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt", sep="/"), header=FALSE, as.is=TRUE)
keep7 <- !(featureNames(MSetsq_28Flt) %in% as.character(multimapProbes$V1))   # Watch out!!, download and include annotation for EPIC array
MSetsq_28Flt <- MSetsq_28Flt[keep7,]

# Visualize the MDS plots
pdf(file=paste(resultsDirectory, "04-MDS_plots_pre_post_filtering_probes_28.pdf", sep = "/"), paper="a4r", height=21, width=28, onefile=TRUE)
par(mfrow=c(1,1))
plotMDS(getM(MSetsq_28), top=1000, gene.selection="common", col=phenoData_28_filtered$color, main="MDS of Noob M-values (Case/Control)", labels=phenoData_28_filtered$Condition)
legend("bottomright", legend=levels(as.factor(phenoData_28_filtered$Sample_Group)), text.col=levels(as.factor(phenoData_28_filtered$color)), bg="white", cex=0.8)
plotMDS(getM(MSetsq_28Flt), top=1000, gene.selection="common", col=phenoData_28_filtered$color, dim=c(1,3), main="MDS Noob post-filtering", labels=phenoData_28_filtered$Condition)
legend("bottomright", legend=levels(as.factor(phenoData_28_filtered$Sample_Group)), text.col=levels(as.factor(phenoData_28_filtered$color)), bg="white", cex=0.8)

dev.off()

# calculate BMIQ for correction of typeI and typeII probes:
bVals_28 <- getBeta(MSetsq_28Flt)
bVals_28 <- champ.norm(beta = bVals_28, method = "BMIQ", plotBMIQ=FALSE, arraytype="EPIC", cores=4)
colnames(bVals_28) <- phenoData_28_filtered$Condition

# Calculate Mvalues for statistical analyses:
mVals_28 <- beta2m(bVals_28)

# Write and save M and B values:
setwd(tablesDirectory)
                                                                                                                                                                                                      
# In our case (nanomaterials project), m values will not be used because we can not make statistics with one sample per group
write.table(mVals_28, file = "mVals_28.txt", sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)
# This steps would be to keep groups_28 for further comparisons
#write.table(mVals_28[,1:6], file = "~/mVals_28.txt", sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)
#write.table(mVals_28[,7:10], file = "~/mVals_28.txt", sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)
write.table(bVals_28, file = "bVals_28.txt", sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)
# This steps would be to keep groups_28 for further comparisons
#write.table(bVals_28[,1:6], file = "~/bVals_28.txt", sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)
#write.table(bVals_28[,7:10], file = "~/bVals_28.txt", sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)

write.table(phenoData_28_filtered, file = "phenoData_28.txt", sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)

# Plot PCAs
#colnames(bVals_28) <- phenoData_28_filtered$Condition
PCA_28 <- prcomp(t(bVals_28), scale. = FALSE)
percentVar <- round(100*PCA_28$sdev^2/sum(PCA_28$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

pdf(file=paste(resultsDirectory, "05-PCA_plots_conditions_28.pdf", sep = "/"), paper="a4r", height=21, width=28, onefile=TRUE)
fviz_pca_ind(PCA_28, repel = T, geom.ind = c("point", "text"), fill.ind = c(as.character(phenoData_28$Sample_Group)), col.ind = "black", pointshape = 21, pointsize = 4, addEllipses = F, palette = "jco")

fviz_pca_ind(PCA_28, repel = T, geom.ind = c("point", "text"), fill.ind = c(as.character(phenoData_28$Subclasificacion_histologica)), col.ind = "black", pointshape = 21, pointsize = 4, addEllipses = F, palette = "jco")

fviz_pca_ind(PCA_28, repel = T, geom.ind = c("point", "text"), fill.ind = c(as.character(phenoData_28$IHQ_FT)), col.ind = "black", pointshape = 21, pointsize = 4, addEllipses = F, palette = "jco")

fviz_pca_ind(PCA_28, repel = T, geom.ind = c("point", "text"), fill.ind = c(as.character(phenoData_28$Fallecido)), col.ind = "black", pointshape = 21, pointsize = 4, addEllipses = F, palette = "jco")

fviz_pca_ind(PCA_28, repel = T, geom.ind = c("point", "text"), fill.ind = c(as.character(phenoData_28$Comentarios)), col.ind = "black", pointshape = 21, pointsize = 4, addEllipses = F, palette = "jco")

fviz_pca_ind(PCA_28, repel = T, geom.ind = c("point", "text"), fill.ind = c(as.character(phenoData_28$Comorb_HTA)), col.ind = "black", pointshape = 21, pointsize = 4, addEllipses = F, palette = "jco")

fviz_pca_ind(PCA_28, repel = T, geom.ind = c("point", "text"), fill.ind = c(as.character(phenoData_28$Comorb_Hipopituitarismo)), col.ind = "black", pointshape = 21, pointsize = 4, addEllipses = F, palette = "jco")

fviz_pca_ind(PCA_28, repel = T, geom.ind = c("point", "text"), fill.ind = c(as.character(phenoData_28$Invasion_Seno_Cavernoso)), col.ind = "black", pointshape = 21, pointsize = 4, addEllipses = F, palette = "jco")

fviz_pca_ind(PCA_28, repel = T, geom.ind = c("point", "text"), fill.ind = c(as.character(phenoData_28$Grado_de_invasion_escala_Knosp)), col.ind = "black", pointshape = 21, pointsize = 4, addEllipses = F, palette = "jco")

fviz_pca_ind(PCA_28, repel = T, geom.ind = c("point", "text"), fill.ind = c(as.character(phenoData_28$Fecha_de_diagnostico)), col.ind = "black", pointshape = 21, pointsize = 4, addEllipses = F, palette = "jco")

dev.off()

# Observamos que la muestra "29_No" pertenece al grupo de los 4 fallecidos y es la única
# con Transplante renal (Comentarios) y con IHQ_FT = NULO. 
# La muestra "11_No" es la única con Subclasificación histológica: CTS.