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

### Annotate CpG sites from EPIC array:

library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Set up the path for the project location:
basedir <- "/media/rauldiul/Expansion/1_ROCIOyJAVI/RG_Macroadenomas"
tablesDirectory <- file.path(basedir, "Tables/")

# Set working directory:
setwd(tablesDirectory)

# Load annotationEPIC and bValues
annotationEPIC <- as.data.frame(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19))
rownames(annotationEPIC) <- annotationEPIC$Name

bVals <- read.csv("bVals_28.txt", sep = "\t", stringsAsFactors = F)

# Filter annotationEPIC according to the CpGs in bVals, customize and save the table
annotationEPIC <- annotationEPIC[rownames(bVals),]
annotationEPIC$chr <- as.character(annotationEPIC$chr)
annotationEPIC$Name <- as.character(annotationEPIC$Name)
annotationEPIC$pos <- as.integer(annotationEPIC$pos)
annotationEPIC$end <- as.integer(annotationEPIC$pos + 1)
annotationEPIC <- annotationEPIC[order(annotationEPIC[,"chr"], annotationEPIC[,"pos"]),]
annotationEPIC <- annotationEPIC[,c("chr", "pos", "end", "strand", "Name", "Type", "Probe_rs", "Probe_maf", "Islands_Name", "Relation_to_Island", "DMR", "Phantom4_Enhancers", "Phantom5_Enhancers", "X450k_Enhancer", "HMM_Island", "GencodeCompV12_NAME", "GencodeCompV12_Accession", "GencodeCompV12_Group", "Methyl27_Loci", "Methyl450_Loci")]
write.table(annotationEPIC, file = "annotation_EPIC.txt", sep = "\t", row.names = F, col.names = T, quote = F)

# Add
library(ChIPseeker)

peak_all <- readPeakFile("annotation_EPIC.txt")
peakAnnoall <- annotatePeak(peak_all, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
all850k <- peakAnnoall@anno@elementMetadata
rownames(all850k) <- as.character(all850k$Name)
all850k$annotation <- gsub("Exon.*", "Exon", all850k$annotation)
all850k$annotation <- gsub("Intron.*", "Intron", all850k$annotation)
all850k$annotation <- gsub("Downstream.*", "Downstream", all850k$annotation)
all850k$annotation <- gsub("Promoter \\(1-2kb\\)", "Distal promoter", all850k$annotation)
all850k$annotation <- gsub("Promoter \\(2-3kb\\)", "Distal promoter", all850k$annotation)

all <- cbind(annotationEPIC[,c(1:3)], all850k)
all <- as.data.frame(all)

# Modify the "Relation_to_Island" column to (later) annotate CpG islands
all$Relation_to_Island <- gsub(pattern = "N_Shelf", replacement = "Shelf", x = all$Relation_to_Island)
all$Relation_to_Island <- gsub(pattern = "S_Shelf", replacement = "Shelf", x = all$Relation_to_Island)
all$Relation_to_Island <- gsub(pattern = "N_Shore", replacement = "Shore", x = all$Relation_to_Island)
all$Relation_to_Island <- gsub(pattern = "S_Shore", replacement = "Shore", x = all$Relation_to_Island)

# Save the table
write.table(all, file = "annotation_EPIC.txt", sep = "\t", col.names = T, row.names = F, quote = F)


###########################   TO ANNOTATE DMPs   ##########################

# Load the annotation_EPIC table
annotation_EPIC <- read.csv("annotation_EPIC.txt", sep = "\t", stringsAsFactors = F)
rownames(annotation_EPIC) <- annotation_EPIC$Name

# For the annotation of the geneId with the gene symbol:
library("biomaRt")
mart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl", host="www.ensembl.org")

# Load the datasets:
R_vs_No <- read.csv("DMPs_R_UnAdjP05_b20.txt", sep = "\t")
Recid_vs_No <- read.csv("DMPs_R1_UnAdjP05_b20.txt", sep = "\t")
RECUR_vs_No <- read.csv("DMPs_R2_UnAdjP05_b20.txt", sep = "\t")
RECUR_parcial <- read.csv("DMPs_R2_parcial_UnAdjP05_b20.txt", sep = "\t")

# Annotate datasets:
# Full annotation of datasets including statistics (p-values...) and Gene Symbol
R_vs_No$Name <- rownames(R_vs_No)
R_vs_No_all <- merge(R_vs_No, annotation_EPIC, by = "Name")
symbol_R_vs_No_all <- ensembldb::select(mart, keys = R_vs_No_all$geneId, keytype = "entrezgene_id", columns = c("entrezgene_id", "external_gene_name"))
R_vs_No_all$symbol <- symbol_R_vs_No_all$external_gene_name[match(R_vs_No_all$geneId, symbol_R_vs_No_all$entrezgene_id)]
R_vs_No_HYPER <- R_vs_No_all[which(R_vs_No_all$Methylstat == "hyper"),]
R_vs_No_hypo <- R_vs_No_all[which(R_vs_No_all$Methylstat == "hypo"),]

Recid_vs_No$Name <- rownames(Recid_vs_No)
Recid_all <- merge(Recid_vs_No, annotation_EPIC, by = "Name")
symbol_Recid_all <- ensembldb::select(mart, keys = Recid_all$geneId, keytype = "entrezgene_id", columns = c("entrezgene_id", "external_gene_name"))
Recid_all$symbol <- symbol_Recid_all$external_gene_name[match(Recid_all$geneId, symbol_Recid_all$entrezgene_id)]
Recid_HYPER <- Recid_all[which(Recid_all$Methylstat == "hyper"),]
Recid_hypo <- Recid_all[which(Recid_all$Methylstat == "hypo"),]

RECUR_vs_No$Name <- rownames(RECUR_vs_No)
RECUR_all <- merge(RECUR_vs_No, annotation_EPIC, by = "Name")
symbol_RECUR_all <- ensembldb::select(mart, keys = RECUR_all$geneId, keytype = "entrezgene_id", columns = c("entrezgene_id", "external_gene_name"))
RECUR_all$symbol <- symbol_RECUR_all$external_gene_name[match(RECUR_all$geneId, symbol_RECUR_all$entrezgene_id)]
RECUR_HYPER <- RECUR_all[which(RECUR_all$Methylstat == "hyper"),]
RECUR_hypo <- RECUR_all[which(RECUR_all$Methylstat == "hypo"),]

RECUR_parcial$Name <- rownames(RECUR_parcial)
RECUR_parcial_all <- merge(RECUR_parcial, annotation_EPIC, by = "Name")
symbol_RECUR_parcial_all <- ensembldb::select(mart, keys = RECUR_parcial_all$geneId, keytype = "entrezgene_id", columns = c("entrezgene_id", "external_gene_name"))
RECUR_parcial_all$symbol <- symbol_RECUR_parcial_all$external_gene_name[match(RECUR_parcial_all$geneId, symbol_RECUR_parcial_all$entrezgene_id)]
RECUR_parcial_HYPER <- RECUR_parcial_all[which(RECUR_parcial_all$Methylstat == "hyper"),]
RECUR_parcial_hypo <- RECUR_parcial_all[which(RECUR_parcial_all$Methylstat == "hypo"),]


#################
# Write the resulting dataframes in txt files:
write.table(R_vs_No_all, file = paste(tablesDirectory, "Annotated_R_vs_No_all.txt", sep = "/"), sep = "\t", row.names = TRUE)
write.table(R_vs_No_HYPER, file = paste(tablesDirectory, "Annotated_R_vs_No_HYPER.txt", sep = "/"), sep = "\t", row.names = TRUE)
write.table(R_vs_No_hypo, file = paste(tablesDirectory, "Annotated_R_vs_No_hypo.txt", sep = "/"), sep = "\t", row.names = TRUE)

write.table(Recid_all, file = paste(tablesDirectory, "Annotated_Recid_all.txt", sep = "/"), sep = "\t", row.names = TRUE)
write.table(Recid_HYPER, file = paste(tablesDirectory, "Annotated_Recid_HYPER.txt", sep = "/"), sep = "\t", row.names = TRUE)
write.table(Recid_hypo, file = paste(tablesDirectory, "Annotated_Recid_hypo.txt", sep = "/"), sep = "\t", row.names = TRUE)

write.table(RECUR_all, file = paste(tablesDirectory, "Annotated_RECUR_all.txt", sep = "/"), sep = "\t", row.names = TRUE)
write.table(RECUR_HYPER, file = paste(tablesDirectory, "Annotated_RECUR_HYPER.txt", sep = "/"), sep = "\t", row.names = TRUE)
write.table(RECUR_hypo, file = paste(tablesDirectory, "Annotated_RECUR_hypo.txt", sep = "/"), sep = "\t", row.names = TRUE)

write.table(RECUR_parcial_all, file = paste(tablesDirectory, "Annotated_RECUR_parcial_all.txt", sep = "/"), sep = "\t", row.names = TRUE)
write.table(RECUR_parcial_HYPER, file = paste(tablesDirectory, "Annotated_RECUR_parcial_HYPER.txt", sep = "/"), sep = "\t", row.names = TRUE)
write.table(RECUR_parcial_hypo, file = paste(tablesDirectory, "Annotated_RECUR_parcial_hypo.txt", sep = "/"), sep = "\t", row.names = TRUE)


# Beta values and M-values
#mVals <- as.data.frame(mVals[annotation_EPIC$Name,])
#write.table(mVals, file = "mVals.txt", sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)

#bVals <- as.data.frame(bVals[annotation_EPIC$Name,])
#write.table(bVals, file = "bVals.txt", sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)
