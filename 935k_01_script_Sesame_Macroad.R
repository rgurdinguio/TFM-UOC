#######################################################################
# PROJECT: Non-functioning Pituitary Macroadenomas
#######################################################################
## Basic Preprocessing with Sesame
## Script for the analysis of Methylation array V2 data (EPICv2 935K)
## Author: Rocío G. Urdinguio
## (adapted from Sesame info: [sesame github vignette](https://github.com/zwdzwd/sesame/blob/devel/vignettes/sesame.Rmd) and [Sesame bioconductor vignette](https://bioconductor.org/packages/release/bioc/vignettes/sesame/inst/doc/sesame.html)).

#######################################################################
### Preparation of Sesame package 
#######################################################################

# Install and library Sesame.
#BiocManager::install("sesame")
library(sesame)
library(SummarizedExperiment)
#library(sesameData)

## As sesame and sesameData are under active development, this documentation is
## specific to the following version of R, sesame, sesameData and ExperimentHub:
sesame_checkVersion()
#sesameDataCache() #Must be run only after a NEW INSTALATION

#######################################################################
### Required packages and directories
#######################################################################

# Set up the path for the project:
getwd()
# If we are in "scr" folder, move up a folder by running:
setwd("..")
getwd()
basedir <- getwd()

#IDATs localization:
idat_dir <- file.path(basedir, "Raw_Data_935K/")

# Set up the directory for the results:
results_dir <- file.path(basedir, "PDF/")
tables_dir <- file.path(basedir, "Tables/")

#Load the required packages:
#library(stats)
# Load the required packages:
#library(sesame)
#library(SummarizedExperiment)
#library(dplyr)
#library(tidyr)
#library(readxl)
#library(factoextra) # for PCA analysis
#library(wheatmap)
#library(ggplot2)
#library(gprofiler2)

#######################################################################
### GET BETAS (without prep or collapse) -- NOT FOR US
#######################################################################

## GET BETAS
betas = openSesame(idat_dir, BPPARAM = BiocParallel::SnowParam(4), func = getBetas) 
# getBetas is the default (and could be omitted)

## RETURN THE `SigDF` LIST:  
sdfs = openSesame(idat_dir, BPPARAM = BiocParallel::SnowParam(4), func = NULL) 
#`SigDF` is a data.frame with 7 columns: Probe_ID, MG, MR, UG, UR, col and mask. 
#The `col` column specifies the color channel and takes G, R and 2: The Infinium-I probes carry G and R in col to indicate the designed color, while the Infinium-II probes carry 2.

## SUMMARIZE RESULTING SigDF:
sesameQC_calcStats(sdfs[[1]], "numProbes") #para la primera muestra

## CALCULATE SNP ALLELE FREQUENCIES:  
allele_freqs = openSesame(idat_dir, BPPARAM = BiocParallel::SnowParam(4), func = getAFs)

## GENERATE THE DETECTION p-values (e.g., for GEO upload):  
pvals = openSesame(idat_dir, BPPARAM = BiocParallel::SnowParam(4), func = pOOBAH,
                   return.pval=TRUE) 

## WRITE AND SAVE TABLES:  
write.table(betas, file = paste(tables_dir, "betas.txt", sep = "/"), sep = "\t",
            col.names = TRUE, row.names = TRUE, quote = F)
sdf_write_table(sdfs, file= paste(tables_dir, "sdfs.tsv", sep = "/"), sep="\t",
                quote=FALSE)
write.table(allele_freqs, file = paste(tables_dir, "allele_freqs.txt", sep = "/"), 
            sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)
write.table(pvals, file = paste(tables_dir, "pvals.txt", sep = "/"), sep = "\t",
            col.names = TRUE, row.names = TRUE, quote = F)

#######################################################################
#### Preprocessing Function Code {#prep}  
#######################################################################

# DATA PREPROCESSING:  
sdfs_prep = openSesame(sdfs, BPPARAM = BiocParallel::SnowParam(4), func=NULL, prep="QCDPB")
betas_prep = openSesame(idat_dir, BPPARAM = BiocParallel::SnowParam(4), func = getBetas,
                        prep="QCDPB") 
allele_freqs_prep = openSesame(idat_dir, BPPARAM = BiocParallel::SnowParam(4), func = getAFs,
                               prep="QCDPB") 
pvals_prep = openSesame(idat_dir, BPPARAM = BiocParallel::SnowParam(4), func = pOOBAH, 
                        return.pval=TRUE, prep="QCDPB") 

# WRITE AND SAVE TABLES:  
write.table(betas_prep, file = paste(tables_dir, "betas_prep.txt", sep = "/"), 
            sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)
sdf_write_table(sdfs_prep, file= paste(tables_dir, "sdfs_prep.tsv", sep = "/"), 
                sep="\t", quote=FALSE)
write.table(allele_freqs_prep, file = paste(tables_dir, "allele_freqs_prep.txt", sep = "/"),
            sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)
write.table(pvals_prep, file = paste(tables_dir, "pvals_prep.txt", sep = "/"), 
            sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)

#######################################################################
### Collapse measurements to cg prefixes  
#######################################################################

# Collapse measurements to cg prefixes 
betas_collapsed = openSesame(idat_dir, BPPARAM = BiocParallel::SnowParam(4), 
                             func = getBetas, collapseToPfx = TRUE, 
                             collapseMethod = "mean")

################### THIS IS THE ONE WE WILL USE: #####################
betas_prep_collapsed = openSesame(idat_dir, BPPARAM = BiocParallel::SnowParam(4), 
                                  func = getBetas, prep="QCDPB", collapseToPfx = TRUE,
                                  collapseMethod = "mean")

#By default the method for collapsing is to make means, but it can also be switched to `min detection p-value`:   
betas_collapsed_minPval = openSesame(idat_dir, BPPARAM = BiocParallel::SnowParam(4), 
                                     func = getBetas, collapseToPfx = TRUE, 
                                     collapseMethod = "minPval")

# WRITE AND SAVE TABLES:   
write.table(betas_collapsed, file = paste(tables_dir, "betas_collapsed.txt", sep = "/"),
            sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)
write.table(betas_prep_collapsed, file = paste(tables_dir, "betas_prep_collapsed.txt", sep = "/"),
            sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)
write.table(betas_collapsed_minPval, file = paste(tables_dir, "betas_collapsed_minPval.txt", sep = "/"),
            sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)

### Save Environment  
save.image(paste(basedir, "scr/01_SCRIPT_SUMMARY_sesame_Macroad935_Environment.RData", sep="/"))

#######################################################################
## QC after Basic Preprocessing with Sesame
#######################################################################

# CALCULATE QUALITY METRICS
# calculate metrics on all IDATs in a specific folder
qcs = openSesame(idat_dir, prep="", func=sesameQC_calcStats)
# COMBINE A LIST OF sesameQC into a data frame:
qcs_df <- do.call(rbind, lapply(qcs, as.data.frame))
# WRITE AND SAVE:
write.table(qcs_df, file = paste(tables_dir,"qcs_df.txt", sep = "/"), sep = "\t", 
            col.names = TRUE, row.names = TRUE, quote = F)

### Signal detection
# We consider signal detection the most important QC metric.
qc_detection = openSesame(sdfs, prep="", func=sesameQC_calcStats, funs="detection")
qc_detection_df <- do.call(rbind, lapply(qc_detection, as.data.frame))

### Other Quality Control Plots
# `sesameQC_plotBar()` takes a list of sesameQC objects and creates bar plot for each metric calculated.
# `sesameQC_plotRedGrnQQ()` graphs the dye bias between the two color channels.
# `sesameQC_plotIntensVsBetas()` plots the relationship between $\beta$ values and signal intensity and can be used to diagnose artificial readout and influence of signal background.
# `sesameQC_plotHeatSNPs()` plots SNP probes and can be used to detect sample swaps.

#### QC Stats Bar plot
# The fraction of detection failures are signs of masking due to variety of reasons including failed detection, high background, putative low quality probes etc. 
pdf(file=paste(results_dir, "QC_1_plotBar.pdf", sep = "/"), paper="a4r", 
    height=21, width=28, onefile=TRUE)
sesameQC_plotBar(lapply(sdfs, sesameQC_calcStats, "detection"))
sesameQC_plotBar(lapply(sdfs_prep, sesameQC_calcStats, "detection"))
dev.off()

#### Dye bias Q-Q plot
# Dye bias is shown by an off-diagonal q-q plot of the red (x-axis) and green signal (y-axis).
pdf(file=paste(results_dir, "QC_2_plotRedGrnQQ.pdf", sep = "/"), paper="a4r", height=21,
     width=28, onefile=TRUE)
for (i in 1:length(sdfs_prep)) {
  sesameQC_plotRedGrnQQ(sdfs[[i]])
  sesameQC_plotRedGrnQQ(sdfs_prep[[i]])} 
dev.off()
# Formatting: 
pdf(file=paste(results_dir, "QC_2_plotRedGrnQQ_v3.pdf", sep = "/"), paper="a4r", height=21,
   width=28) #, onefile=TRUE)
par(mfrow = c(5, 7), cex=0.75, mar = c(3.1, 2.1, 2.1, 1.1)) #lwd=0.5,
for (i in 1:length(sdfs)) {
 sesameQC_plotRedGrnQQ(sdfs[[i]], cex.lab=0.75, cex.axis=0.75, cex.main=0.75)} 
par(mfrow = c(5, 7))
for (i in 1:length(sdfs_prep)) {
 sesameQC_plotRedGrnQQ(sdfs_prep[[i]], cex.lab=0.75, cex.axis=0.75, main=0.75)} 
dev.off()

#### Intensity-beta plot  
pdf(file=paste(results_dir, "QC_3_plotIntensVsBetas.pdf", sep = "/"), paper="a4r", height=21, width=28, onefile=TRUE)
for (i in 1:length(sdfs_prep)) {
 sesameQC_plotIntensVsBetas(sdfs[[i]])
 sesameQC_plotIntensVsBetas(sdfs_prep[[i]])}
dev.off()

# Otra disposición
pdf(file=paste(results_dir, "QC_3_plotIntensVsBetas_v2.pdf", sep = "/"), paper="a4r", height=21, width=28) #, onefile=TRUE)
par(mfrow = c(2, 2), cex = 0.9)
for (i in 1:length(sdfs_prep)) {
 sesameQC_plotIntensVsBetas(sdfs[[i]], cex.lab=0.9, cex.axis=0.9, cex.main=0.9)
 sesameQC_plotIntensVsBetas(sdfs_prep[[i]], cex.lab=0.9, cex.axis=0.9, cex.main=0.9)}
dev.off()
# Beta value is more influenced by signal background for probes with low signal intensities. 
# The following plot shows this dependency and the extent of probes with low signal intensity.
   
#### Genotype validation  
# Extra SNP allele frequencies can be obtained by summing up methylated and umethylated alleles of color-channel-switching probes.   
# These allele frequencies can be combined with explicit SNP probes:  
head(getAFs(sdfs[[1]])) 
#Extract explicit and Infinium-I-derived SNPs to identify potential sample swaps:
pdf(file=paste(results_dir, "QC_4_plotHeatSNPs.pdf", sep = "/"), paper="a4r", height=21, width=28, onefile=TRUE)
sesameQC_plotHeatSNPs(sdfs)
sesameQC_plotHeatSNPs(sdfs_prep)
dev.off()
# One can also output the allele frequencies and output a VCF file with genotypes. 
# This requires additional SNP information (ref and alt alleles), which can be downloaded using the following code:
#head(formatVCF(sdf)) #Atención: Habria que preparalo para que funcione con el EPICv2
#One can output to actual VCF file with a header by formatVCF(sdf, vcf=path_to_vcf) #Atención: Habria que preparalo para que funcione con el EPICv2

#### Bisulfite conversion
#bisConversionControl(sdfs[[1]]) #Atencion: Habria que preparalo para que funcione con el EPICv2

### Save Environment
save.image(paste(basedir, "scr/01_SCRIPT_SUMMARY_sesame_Macroad935_Environment.RData", sep="/"))
