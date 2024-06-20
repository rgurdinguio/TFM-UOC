#######################################################################
# PROJECT: Non-functioning Pituitary Macroadenomas
#######################################################################
### Script for the analysis of Methylation array V2 data (EPICv2 935K)
### 
### Author: Rocío G. Urdinguio
########################################################################

###############################################################
###      Generates barplots for loci representation         ###
###############################################################

library(RColorBrewer)
library(ggplot2)
library(tidyr)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(minfi)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene


#################################################################
### SET UP THE PATHS AND DIRECTORIES
#################################################################

# Set up basedir
getwd()
setwd("..")
getwd()
basedir <- getwd()

# Set up the directory for the results:
results_dir <- file.path(basedir, "PDF/")
tables_dir <- file.path(basedir, "Tables/")

##################################################################
### LOAD BETAS AND phenodata
##################################################################

# Load betas
betas <- fread(file.path(tables_dir, "02_betas_noob.txt.gz"))
betas.flt <- fread(file.path(tables_dir, "02_betas_noob_filtered.txt.gz"))
# Also array annotation
annotationEPIC <- fread(file.path(tables_dir, "02_EPICv2_annotation.txt.gz"))
#annotationEPIC <- read.csv("annotationEPIC.txt", sep = "\t", stringsAsFactors = F)
#rownames(annotationEPIC) <- annotationEPIC$ID2
#Rocío
phenodata <- fread(file.path(tables_dir,"02_phenodata.txt"))

###########################   TO ANNOTATE DMPs   ##########################

##################################################################
### LOAD Datasets
##################################################################

# Load the annotation_EPIC table
annotationEPIC <- fread(file.path(tables_dir, "02_EPICv2_annotation.txt.gz"))
#rownames(annotation_EPIC) <- annotation_EPIC$ID2

# For the annotation of the geneId with the gene symbol:
library("biomaRt")
mart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl", host="https://www.ensembl.org")

# Load the datasets:
RvsNo_UnAdjP005_b10 <- read.csv(paste(tables_dir, "test_result_RvsNo_UnAdjP005_b10.txt", sep = "/"), 
                             sep = "\t") #, stringsAsFactors = F)
# Annotate datasets:
# Full annotation of datasets including statistics (p-values...) and Gene Symbol
#BiocManager::install("ensembldb")
library(ensembldb)
RvsNo_UnAdjP005_b10$ID2 <- RvsNo_UnAdjP005_b10$Probe_ID
R_vs_No_all <- merge(RvsNo_UnAdjP005_b10, annotationEPIC, by = "ID2")
symbol_R_vs_No_all <- ensembldb::select(mart, keys = R_vs_No_all$chipseeker_geneId, keytype = "entrezgene_id", columns = c("entrezgene_id", "external_gene_name"))
R_vs_No_all$symbol <- symbol_R_vs_No_all$external_gene_name[match(R_vs_No_all$chipseeker_geneId, symbol_R_vs_No_all$entrezgene_id)]
R_vs_No_HYPER <- R_vs_No_all[which(R_vs_No_all$Methylstat == "hyper"),]
R_vs_No_hypo <- R_vs_No_all[which(R_vs_No_all$Methylstat == "hypo"),]

#################
# Write the resulting dataframes in txt files:
write.table(R_vs_No_all, file = paste(tables_dir, "Annotated_R_vs_No_all.txt", sep = "/"), sep = "\t", row.names = TRUE)
write.table(R_vs_No_HYPER, file = paste(tables_dir, "Annotated_R_vs_No_HYPER.txt", sep = "/"), sep = "\t", row.names = TRUE)
write.table(R_vs_No_hypo, file = paste(tables_dir, "Annotated_R_vs_No_hypo.txt", sep = "/"), sep = "\t", row.names = TRUE)

###################################  Number of CpGs  ########################################

numbers_CpGs_Back <- c(length(annotationEPIC$ID2), 0, 0)

numbers_CpGs_R <- c((length(R_vs_No_HYPER$ID2) + length(R_vs_No_hypo$ID2)),
                     length(R_vs_No_HYPER$ID2),
                     length(R_vs_No_hypo$ID2))

####  OR  ####
numbers_CpGs_Total <- c(length(annotationEPIC$ID2), 
                        (length(R_vs_No_HYPER$ID2) + length(R_vs_No_hypo$ID2)))

numbers_CpGs_HYPER <- c(0,
                        length(R_vs_No_HYPER$ID2))

numbers_CpGs_Hypo <- c(0,
                       length(R_vs_No_hypo$ID2))

Comparison <- c("Back", "R_vs_No")
Methylstat <- c("Total", "Hyper","Hypo")

df_CpG_counts <- data.frame(numbers_CpGs_Total, numbers_CpGs_HYPER, numbers_CpGs_Hypo)
colnames(df_CpG_counts) <- Methylstat
row.names(df_CpG_counts) <- Comparison

# Write the resulting dataframe in a CSV file:
write.table(df_CpG_counts, file = paste(tables_dir, "CpG_counts_UnAdjP005_b10.csv", sep = "/"), sep = ",", row.names = TRUE)

#######################################  ISLANDS  ############################################

# Establish the data frames and perform the plots
group_island_R <- c(rep("Back",4), rep("HYPER", 4), rep("hypo", 4))
class_island <- c(rep(c("CGI;OpenSea", "CGI;Shore", "CGI;Shelf", "CGI;Island"), 3))

numbers_island_R <- c(length(annotationEPIC$Relation_to_CpG_Island[which(annotationEPIC$Relation_to_CpG_Island == "CGI;OpenSea")]), length(annotationEPIC$Relation_to_CpG_Island[which(annotationEPIC$Relation_to_CpG_Island == "CGI;Shore")]), length(annotationEPIC$Relation_to_CpG_Island[which(annotationEPIC$Relation_to_CpG_Island == "CGI;Shelf")]), length(annotationEPIC$Relation_to_CpG_Island[which(annotationEPIC$Relation_to_CpG_Island == "CGI;Island")]),
                      length(R_vs_No_HYPER$Relation_to_CpG_Island[which(R_vs_No_HYPER$Relation_to_CpG_Island == "CGI;OpenSea")]), length(R_vs_No_HYPER$Relation_to_CpG_Island[which(R_vs_No_HYPER$Relation_to_CpG_Island == "CGI;Shore")]), length(R_vs_No_HYPER$Relation_to_CpG_Island[which(R_vs_No_HYPER$Relation_to_CpG_Island == "CGI;Shelf")]), length(R_vs_No_HYPER$Relation_to_CpG_Island[which(R_vs_No_HYPER$Relation_to_CpG_Island == "CGI;Island")]),
                      length(R_vs_No_hypo$Relation_to_CpG_Island[which(R_vs_No_hypo$Relation_to_CpG_Island == "CGI;OpenSea")]), length(R_vs_No_hypo$Relation_to_CpG_Island[which(R_vs_No_hypo$Relation_to_CpG_Island == "CGI;Shore")]), length(R_vs_No_hypo$Relation_to_CpG_Island[which(R_vs_No_hypo$Relation_to_CpG_Island == "CGI;Shelf")]), length(R_vs_No_hypo$Relation_to_CpG_Island[which(R_vs_No_hypo$Relation_to_CpG_Island == "CGI;Island")]))

df_island_R <- data.frame(group_island_R,class_island,numbers_island_R)
df_island_R$class_island <- factor(df_island_R$class_island, levels = c("CGI;Island", "CGI;Shelf", "CGI;Shore", "CGI;OpenSea"))
df_island_R$group_island_R <- factor(df_island_R$group_island_R, levels = c("Back", "HYPER", "hypo"))


pdf(file = paste(results_dir, "DMPs_genomic_distribution_CpGs_v2.pdf", sep = "/"), height = 5, width = 4, onefile = T)
ggplot(df_island_R, aes(x = group_island_R)) + geom_bar(aes(weight= numbers_island_R, fill = class_island), position = 'fill') + scale_fill_manual(values = rev(c("#d09c9c", "#bf7272", "#a72b2b", "#761212"))) + theme_bw() + ggtitle(label = "R_vs_No \nCpG distribution CpGI status") + ylab(label = "Relative Frequency") + theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5)) + xlab("")

#######################################  REGIONS  ############################################

group_region_R <- c(rep("Back",8), rep("HYPER", 8), rep("hypo", 8))
class_region <- c(rep(c("Distal promoter", "Promoter (<=1kb)", "5' UTR", "Intron", "Exon", "3' UTR", "Downstream", "Distal Intergenic"), 3))

numbers_region_R <- c(length(annotationEPIC$chipseeker_annotation[which(annotationEPIC$chipseeker_annotation == "Distal promoter")]), length(annotationEPIC$chipseeker_annotation[which(annotationEPIC$chipseeker_annotation == "Promoter (<=1kb)")]), length(annotationEPIC$chipseeker_annotation[which(annotationEPIC$chipseeker_annotation == "5' UTR")]), length(annotationEPIC$chipseeker_annotation[which(annotationEPIC$chipseeker_annotation == "Intron")]), length(annotationEPIC$chipseeker_annotation[which(annotationEPIC$chipseeker_annotation == "Exon")]), length(annotationEPIC$chipseeker_annotation[which(annotationEPIC$chipseeker_annotation == "3' UTR")]), length(annotationEPIC$chipseeker_annotation[which(annotationEPIC$chipseeker_annotation == "Downstream")]), length(annotationEPIC$chipseeker_annotation[which(annotationEPIC$chipseeker_annotation == "Distal Intergenic")]),
                      length(R_vs_No_HYPER$chipseeker_annotation[which(R_vs_No_HYPER$chipseeker_annotation == "Distal promoter")]), length(R_vs_No_HYPER$chipseeker_annotation[which(R_vs_No_HYPER$chipseeker_annotation == "Promoter (<=1kb)")]), length(R_vs_No_HYPER$chipseeker_annotation[which(R_vs_No_HYPER$chipseeker_annotation == "5' UTR")]), length(R_vs_No_HYPER$chipseeker_annotation[which(R_vs_No_HYPER$chipseeker_annotation == "Intron")]), length(R_vs_No_HYPER$chipseeker_annotation[which(R_vs_No_HYPER$chipseeker_annotation == "Exon")]), length(R_vs_No_HYPER$chipseeker_annotation[which(R_vs_No_HYPER$chipseeker_annotation == "3' UTR")]), length(R_vs_No_HYPER$chipseeker_annotation[which(R_vs_No_HYPER$chipseeker_annotation == "Downstream")]), length(R_vs_No_HYPER$chipseeker_annotation[which(R_vs_No_HYPER$chipseeker_annotation == "Distal Intergenic")]),
                      length(R_vs_No_hypo$chipseeker_annotation[which(R_vs_No_hypo$chipseeker_annotation == "Distal promoter")]), length(R_vs_No_hypo$chipseeker_annotation[which(R_vs_No_hypo$chipseeker_annotation == "Promoter (<=1kb)")]), length(R_vs_No_hypo$chipseeker_annotation[which(R_vs_No_hypo$chipseeker_annotation == "5' UTR")]), length(R_vs_No_hypo$chipseeker_annotation[which(R_vs_No_hypo$chipseeker_annotation == "Intron")]), length(R_vs_No_hypo$chipseeker_annotation[which(R_vs_No_hypo$chipseeker_annotation == "Exon")]), length(R_vs_No_hypo$chipseeker_annotation[which(R_vs_No_hypo$chipseeker_annotation == "3' UTR")]), length(R_vs_No_hypo$chipseeker_annotation[which(R_vs_No_hypo$chipseeker_annotation == "Downstream")]), length(R_vs_No_hypo$chipseeker_annotation[which(R_vs_No_hypo$chipseeker_annotation == "Distal Intergenic")]))

df_region_R <- data.frame(group_region_R,class_region,numbers_region_R)
df_region_R$class_region <- factor(df_region_R$class_region, levels = c("Distal promoter", "Promoter (<=1kb)", "5' UTR", "Intron", "Exon", "3' UTR", "Downstream", "Distal Intergenic"))
df_region_R$group_region_R <- factor(df_region_R$group_region_R, levels = c("Back", "HYPER", "hypo"))

ggplot(df_region_R, aes(x = group_region_R)) + geom_bar(aes(weight= numbers_region_R, fill = class_region), position = 'fill') + scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")) + theme_bw() + ggtitle(label = "R_vs_No \nCpG distribution Region status") + ylab(label = "Relative Frequency") + theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5)) + xlab("")

dev.off()

######################################################################################################################
###########################################           OTHER GRAPHS          ##########################################
######################################################################################################################

# Establish the data frames and perform the plots
group_island_hyper <- c(rep("Back",4), rep("R_vs_No", 4), rep("Recid", 4), rep("RECUR", 4), rep("RECUR_parcial", 4))
group_island_hypo <- c(rep("Back",4), rep("R_vs_No", 4), rep("Recid", 4), rep("RECUR", 4), rep("RECUR_parcial", 4))
class_island_hyper <- c(rep(c("CGI;OpenSea", "CGI;Shore", "CGI;Shelf", "CGI;Island"), 5))
class_island_hypo <- c(rep(c("CGI;OpenSea", "CGI;Shore", "CGI;Shelf", "CGI;Island"), 5))
numbers_island_hyper <- c(length(annotationEPIC$Relation_to_CpG_Island[which(annotationEPIC$Relation_to_CpG_Island == "CGI;OpenSea")]), length(annotationEPIC$Relation_to_CpG_Island[which(annotationEPIC$Relation_to_CpG_Island == "CGI;Shore")]), length(annotationEPIC$Relation_to_CpG_Island[which(annotationEPIC$Relation_to_CpG_Island == "CGI;Shelf")]), length(annotationEPIC$Relation_to_CpG_Island[which(annotationEPIC$Relation_to_CpG_Island == "CGI;Island")]),
                          length(R_vs_No_HYPER$Relation_to_CpG_Island[which(R_vs_No_HYPER$Relation_to_CpG_Island == "CGI;OpenSea")]), length(R_vs_No_HYPER$Relation_to_CpG_Island[which(R_vs_No_HYPER$Relation_to_CpG_Island == "CGI;Shore")]), length(R_vs_No_HYPER$Relation_to_CpG_Island[which(R_vs_No_HYPER$Relation_to_CpG_Island == "CGI;Shelf")]), length(R_vs_No_HYPER$Relation_to_CpG_Island[which(R_vs_No_HYPER$Relation_to_CpG_Island == "CGI;Island")]),
                          length(Recid_HYPER$Relation_to_CpG_Island[which(Recid_HYPER$Relation_to_CpG_Island == "CGI;OpenSea")]), length(Recid_HYPER$Relation_to_CpG_Island[which(Recid_HYPER$Relation_to_CpG_Island == "CGI;Shore")]), length(Recid_HYPER$Relation_to_CpG_Island[which(Recid_HYPER$Relation_to_CpG_Island == "CGI;Shelf")]), length(Recid_HYPER$Relation_to_CpG_Island[which(Recid_HYPER$Relation_to_CpG_Island == "CGI;Island")]),
                          length(RECUR_HYPER$Relation_to_CpG_Island[which(RECUR_HYPER$Relation_to_CpG_Island == "CGI;OpenSea")]), length(RECUR_HYPER$Relation_to_CpG_Island[which(RECUR_HYPER$Relation_to_CpG_Island == "CGI;Shore")]), length(RECUR_HYPER$Relation_to_CpG_Island[which(RECUR_HYPER$Relation_to_CpG_Island == "CGI;Shelf")]), length(RECUR_HYPER$Relation_to_CpG_Island[which(RECUR_HYPER$Relation_to_CpG_Island == "CGI;Island")]),
                          length(RECUR_parcial_HYPER$Relation_to_CpG_Island[which(RECUR_parcial_HYPER$Relation_to_CpG_Island == "CGI;OpenSea")]), length(RECUR_parcial_HYPER$Relation_to_CpG_Island[which(RECUR_parcial_HYPER$Relation_to_CpG_Island == "CGI;Shore")]), length(RECUR_parcial_HYPER$Relation_to_CpG_Island[which(RECUR_parcial_HYPER$Relation_to_CpG_Island == "CGI;Shelf")]), length(RECUR_parcial_HYPER$Relation_to_CpG_Island[which(RECUR_parcial_HYPER$Relation_to_CpG_Island == "CGI;Island")]))

numbers_island_hypo <- c(length(annotationEPIC$Relation_to_CpG_Island[which(annotationEPIC$Relation_to_CpG_Island == "CGI;OpenSea")]), length(annotationEPIC$Relation_to_CpG_Island[which(annotationEPIC$Relation_to_CpG_Island == "CGI;Shore")]), length(annotationEPIC$Relation_to_CpG_Island[which(annotationEPIC$Relation_to_CpG_Island == "CGI;Shelf")]), length(annotationEPIC$Relation_to_CpG_Island[which(annotationEPIC$Relation_to_CpG_Island == "CGI;Island")]),
                         length(R_vs_No_Hypo$Relation_to_CpG_Island[which(R_vs_No_Hypo$Relation_to_CpG_Island == "CGI;OpenSea")]), length(R_vs_No_Hypo$Relation_to_CpG_Island[which(R_vs_No_Hypo$Relation_to_CpG_Island == "CGI;Shore")]), length(R_vs_No_Hypo$Relation_to_CpG_Island[which(R_vs_No_Hypo$Relation_to_CpG_Island == "CGI;Shelf")]), length(R_vs_No_Hypo$Relation_to_CpG_Island[which(R_vs_No_Hypo$Relation_to_CpG_Island == "CGI;Island")]),
                         length(Recid_Hypo$Relation_to_CpG_Island[which(Recid_Hypo$Relation_to_CpG_Island == "CGI;OpenSea")]), length(Recid_Hypo$Relation_to_CpG_Island[which(Recid_Hypo$Relation_to_CpG_Island == "CGI;Shore")]), length(Recid_Hypo$Relation_to_CpG_Island[which(Recid_Hypo$Relation_to_CpG_Island == "CGI;Shelf")]), length(Recid_Hypo$Relation_to_CpG_Island[which(Recid_Hypo$Relation_to_CpG_Island == "CGI;Island")]),
                         length(RECUR_Hypo$Relation_to_CpG_Island[which(RECUR_Hypo$Relation_to_CpG_Island == "CGI;OpenSea")]), length(RECUR_Hypo$Relation_to_CpG_Island[which(RECUR_Hypo$Relation_to_CpG_Island == "CGI;Shore")]), length(RECUR_Hypo$Relation_to_CpG_Island[which(RECUR_Hypo$Relation_to_CpG_Island == "CGI;Shelf")]), length(RECUR_Hypo$Relation_to_CpG_Island[which(RECUR_Hypo$Relation_to_CpG_Island == "CGI;Island")]),
                         length(RECUR_parcial_Hypo$Relation_to_CpG_Island[which(RECUR_parcial_Hypo$Relation_to_CpG_Island == "CGI;OpenSea")]), length(RECUR_parcial_Hypo$Relation_to_CpG_Island[which(RECUR_parcial_Hypo$Relation_to_CpG_Island == "CGI;Shore")]), length(RECUR_parcial_Hypo$Relation_to_CpG_Island[which(RECUR_parcial_Hypo$Relation_to_CpG_Island == "CGI;Shelf")]), length(RECUR_parcial_Hypo$Relation_to_CpG_Island[which(RECUR_parcial_Hypo$Relation_to_CpG_Island == "CGI;Island")]))

df_island_hyper <- data.frame(group_island_hyper,class_island_hyper,numbers_island_hyper)
df_island_hyper$class_island_hyper <- factor(df_island_hyper$class_island_hyper, levels = c("CGI;Island", "CGI;Shelf", "CGI;Shore", "CGI;OpenSea"))
df_island_hyper$group_island_hyper <- factor(df_island_hyper$group_island_hyper, levels = c("Back", "R_vs_No", "Recid", "RECUR", "RECUR_parcial"))

df_island_hypo <- data.frame(group_island_hypo,class_island_hypo,numbers_island_hypo)
df_island_hypo$class_island_hypo <- factor(df_island_hypo$class_island_hypo, levels = c("CGI;Island", "CGI;Shelf", "CGI;Shore", "CGI;OpenSea"))
df_island_hypo$group_island_hypo <- factor(df_island_hypo$group_island_hypo, levels = c("Back", "R_vs_No", "Recid", "RECUR", "RECUR_parcial"))

setwd(resultsDirectory)

pdf(file = "DMPs_genomic_distribution_CpGs.pdf", height = 5, width = 4, onefile = T)
ggplot(df_island_hyper, aes(x = group_island_hyper)) + geom_bar(aes(weight= numbers_island_hyper, fill = class_island_hyper), position = 'fill') + scale_fill_manual(values = rev(c("#d09c9c", "#bf7272", "#a72b2b", "#761212"))) + theme_bw() + ggtitle(label = "CpG Hyper distribution CpGI status") + ylab(label = "Relative Frequency") + theme(axis.text.x = element_text(angle = 90))
ggplot(df_island_hypo, aes(x = group_island_hypo)) + geom_bar(aes(weight= numbers_island_hypo, fill = class_island_hypo), position = 'fill') + scale_fill_manual(values = rev(c("#d09c9c", "#bf7272", "#a72b2b", "#761212"))) + theme_bw() + ggtitle(label = "CpG Hypo distribution CpGI status") + ylab(label = "Relative Frequency") + theme(axis.text.x = element_text(angle = 90))

###################################################################################

group_region_hyper <- c(rep("Back", 8), rep("R_vs_No", 8))
group_region_hypo <- c(rep("Back", 8), rep("R_vs_No", 8))
class_region_hyper <- c(rep(c("Distal promoter", "Promoter (<=1kb)", "5' UTR", "Intron", "Exon", "3' UTR", "Downstream", "Distal Intergenic"), 2))
class_region_hypo <- c(rep(c("Distal promoter", "Promoter (<=1kb)", "5' UTR", "Intron", "Exon", "3' UTR", "Downstream", "Distal Intergenic"), 2))

numbers_region_hyper <- c(length(annotationEPIC$chipseeker_annotation[which(annotationEPIC$chipseeker_annotation == "Distal promoter")]), length(annotationEPIC$chipseeker_annotation[which(annotationEPIC$chipseeker_annotation == "Promoter (<=1kb)")]), length(annotationEPIC$chipseeker_annotation[which(annotationEPIC$chipseeker_annotation == "5' UTR")]), length(annotationEPIC$chipseeker_annotation[which(annotationEPIC$chipseeker_annotation == "Intron")]), length(annotationEPIC$chipseeker_annotation[which(annotationEPIC$chipseeker_annotation == "Exon")]), length(annotationEPIC$chipseeker_annotation[which(annotationEPIC$chipseeker_annotation == "3' UTR")]), length(annotationEPIC$chipseeker_annotation[which(annotationEPIC$chipseeker_annotation == "Downstream")]), length(annotationEPIC$chipseeker_annotation[which(annotationEPIC$chipseeker_annotation == "Distal Intergenic")]),
                          length(R_vs_No_HYPER$chipseeker_annotation[which(R_vs_No_HYPER$chipseeker_annotation == "Distal promoter")]), length(R_vs_No_HYPER$chipseeker_annotation[which(R_vs_No_HYPER$chipseeker_annotation == "Promoter (<=1kb)")]), length(R_vs_No_HYPER$chipseeker_annotation[which(R_vs_No_HYPER$chipseeker_annotation == "5' UTR")]), length(R_vs_No_HYPER$chipseeker_annotation[which(R_vs_No_HYPER$chipseeker_annotation == "Intron")]), length(R_vs_No_HYPER$chipseeker_annotation[which(R_vs_No_HYPER$chipseeker_annotation == "Exon")]), length(R_vs_No_HYPER$chipseeker_annotation[which(R_vs_No_HYPER$chipseeker_annotation == "3' UTR")]), length(R_vs_No_HYPER$chipseeker_annotation[which(R_vs_No_HYPER$chipseeker_annotation == "Downstream")]), length(R_vs_No_HYPER$chipseeker_annotation[which(R_vs_No_HYPER$chipseeker_annotation == "Distal Intergenic")]))

numbers_region_hypo <- c(length(annotationEPIC$chipseeker_annotation[which(annotationEPIC$chipseeker_annotation == "Distal promoter")]), length(annotationEPIC$chipseeker_annotation[which(annotationEPIC$chipseeker_annotation == "Promoter (<=1kb)")]), length(annotationEPIC$chipseeker_annotation[which(annotationEPIC$chipseeker_annotation == "5' UTR")]), length(annotationEPIC$chipseeker_annotation[which(annotationEPIC$chipseeker_annotation == "Intron")]), length(annotationEPIC$chipseeker_annotation[which(annotationEPIC$chipseeker_annotation == "Exon")]), length(annotationEPIC$chipseeker_annotation[which(annotationEPIC$chipseeker_annotation == "3' UTR")]), length(annotationEPIC$chipseeker_annotation[which(annotationEPIC$chipseeker_annotation == "Downstream")]), length(annotationEPIC$chipseeker_annotation[which(annotationEPIC$chipseeker_annotation == "Distal Intergenic")]),
                         length(R_vs_No_hypo$chipseeker_annotation[which(R_vs_No_hypo$chipseeker_annotation == "Distal promoter")]), length(R_vs_No_hypo$chipseeker_annotation[which(R_vs_No_hypo$chipseeker_annotation == "Promoter (<=1kb)")]), length(R_vs_No_hypo$chipseeker_annotation[which(R_vs_No_hypo$chipseeker_annotation == "5' UTR")]), length(R_vs_No_hypo$chipseeker_annotation[which(R_vs_No_hypo$chipseeker_annotation == "Intron")]), length(R_vs_No_hypo$chipseeker_annotation[which(R_vs_No_hypo$chipseeker_annotation == "Exon")]), length(R_vs_No_hypo$chipseeker_annotation[which(R_vs_No_hypo$chipseeker_annotation == "3' UTR")]), length(R_vs_No_hypo$chipseeker_annotation[which(R_vs_No_hypo$chipseeker_annotation == "Downstream")]), length(R_vs_No_hypo$chipseeker_annotation[which(R_vs_No_hypo$chipseeker_annotation == "Distal Intergenic")]))
                      
df_region_hyper <- data.frame(group_region_hyper,class_region_hyper,numbers_region_hyper)
df_region_hyper$class_region_hyper <- factor(df_region_hyper$class_region_hyper, levels = c("Distal promoter", "Promoter (<=1kb)", "5' UTR", "Intron", "Exon", "3' UTR", "Downstream", "Distal Intergenic")) 
df_region_hyper$group_region_hyper <- factor(df_region_hyper$group_region_hyper, levels = c("Back", "R_vs_No", "Recid", "RECUR", "RECUR_parcial"))

df_region_hypo <- data.frame(group_region_hypo,class_region_hypo,numbers_region_hypo)
df_region_hypo$class_region_hypo <- factor(df_region_hypo$class_region_hypo, levels = c("Distal promoter", "Promoter (<=1kb)", "5' UTR", "Intron", "Exon", "3' UTR", "Downstream", "Distal Intergenic"))  
df_region_hypo$group_region_hypo <- factor(df_region_hypo$group_region_hypo, levels = c("Back", "R_vs_No", "Recid", "RECUR", "RECUR_parcial"))

ggplot(df_region_hyper, aes(x = group_region_hyper)) + geom_bar(aes(weight = numbers_region_hyper, fill = class_region_hyper), position = 'fill') + scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")) + theme_bw() + ggtitle(label = "CpG Hyper distribution Region status") + ylab(label = "Relative Frequency") + theme(axis.text.x = element_text(angle = 90))
ggplot(df_region_hypo, aes(x = group_region_hypo)) + geom_bar(aes(weight = numbers_region_hypo, fill = class_region_hypo), position = 'fill') + scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")) + theme_bw() + ggtitle(label = "CpG Hypo distribution Region status") + ylab(label = "Relative Frequency") + theme(axis.text.x = element_text(angle = 90))

dev.off()
