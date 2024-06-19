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

# Set up the path for the project location:
basedir <- "/media/rauldiul/Expansion/1_ROCIOyJAVI/RG_Macroadenomas"
tablesDirectory <- file.path(basedir, "Tables/")
resultsDirectory <- file.path(basedir, "PDF/")

# Set working directory:
setwd(tablesDirectory)

# Load the datasets:
R_vs_No_all<- read.csv("Annotated_R_vs_No_all.txt", sep = "\t")
R_vs_No_HYPER <- read.csv("Annotated_R_vs_No_HYPER.txt", sep = "\t")
R_vs_No_hypo <- read.csv("Annotated_R_vs_No_hypo.txt", sep = "\t")

Recid_all <- read.csv("Annotated_Recid_all.txt", sep = "\t")
Recid_HYPER <- read.csv("Annotated_Recid_HYPER.txt", sep = "\t")
Recid_hypo <- read.csv("Annotated_Recid_hypo.txt", sep = "\t")

RECUR_all <- read.csv("Annotated_RECUR_all.txt", sep = "\t")
RECUR_HYPER <- read.csv("Annotated_RECUR_HYPER.txt", sep = "\t")
RECUR_hypo <- read.csv("Annotated_RECUR_hypo.txt", sep = "\t")

RECUR_parcial_all <- read.csv("Annotated_RECUR_parcial_all.txt", sep = "\t")
RECUR_parcial_HYPER <- read.csv("Annotated_RECUR_parcial_HYPER.txt", sep = "\t")
RECUR_parcial_hypo <- read.csv("Annotated_RECUR_parcial_hypo.txt", sep = "\t")

annotation_EPIC <- read.csv("annotation_EPIC.txt", sep = "\t", stringsAsFactors = F)
rownames(annotation_EPIC) <- annotation_EPIC$Name

###################################  Number of CpGs  ########################################

numbers_CpGs_Back <- c(length(annotation_EPIC$Name), 0, 0)

numbers_CpGs_R <- c((length(R_vs_No_HYPER$Name) + length(R_vs_No_hypo$Name)),
                     length(R_vs_No_HYPER$Name),
                     length(R_vs_No_hypo$Name))

numbers_CpGs_Recid <- c((length(Recid_HYPER$Name) + length(Recid_hypo$Name)),
                         length(Recid_HYPER$Name),
                         length(Recid_hypo$Name))

numbers_CpGs_RECUR <- c((length(RECUR_HYPER$Name) + length(RECUR_hypo$Name)),
                        length(RECUR_HYPER$Name),
                        length(RECUR_hypo$Name))
                    
numbers_CpGs_RECUR_parcial <- c((length(RECUR_parcial_HYPER$Name) + length(RECUR_parcial_hypo$Name)),
                                 length(RECUR_parcial_HYPER$Name),
                                 length(RECUR_parcial_hypo$Name))

####  OR  ####


numbers_CpGs_Total <- c(length(annotation_EPIC$Name), 
                        (length(R_vs_No_HYPER$Name) + length(R_vs_No_hypo$Name)),
                        (length(Recid_HYPER$Name) + length(Recid_hypo$Name)),
                        (length(RECUR_HYPER$Name) + length(RECUR_hypo$Name)),
                        (length(RECUR_parcial_HYPER$Name) + length(RECUR_parcial_hypo$Name)))

numbers_CpGs_HYPER <- c(0,
                        length(R_vs_No_HYPER$Name),
                        length(Recid_HYPER$Name),
                        length(RECUR_HYPER$Name),
                        length(RECUR_parcial_HYPER$Name))

numbers_CpGs_Hypo <- c(0,
                       length(R_vs_No_hypo$Name),
                       length(Recid_hypo$Name),
                       length(RECUR_hypo$Name),
                       length(RECUR_parcial_hypo$Name))

Comparison <- c("Back", "R_vs_No", "Recid_vs_No", "RECURR_vs_No", "RECURR_vs_Parcial")
Methylstat <- c("Total", "Hyper","Hypo")

df_CpG_counts <- data.frame(numbers_CpGs_Total, numbers_CpGs_HYPER, numbers_CpGs_Hypo)
colnames(df_CpG_counts) <- Methylstat
row.names(df_CpG_counts) <- Comparison

# Write the resulting dataframe in a CSV file:
write.table(df_CpG_counts, file = paste(tablesDirectory, "CpG_counts_UnAdjP05_b20.csv", sep = "/"), sep = ",", row.names = TRUE)


#######################################  ISLANDS  ############################################

# Establish the data frames and perform the plots
group_island_R <- c(rep("Back",4), rep("HYPER", 4), rep("hypo", 4))
group_island_Recid <- c(rep("Back",4), rep("HYPER", 4), rep("hypo", 4))
group_island_RECUR <- c(rep("Back",4), rep("HYPER", 4), rep("hypo", 4))
group_island_RECUR_parcial <- c(rep("Back",4), rep("HYPER", 4), rep("hypo", 4))
class_island <- c(rep(c("OpenSea", "Shore", "Shelf", "Island"), 3))

numbers_island_R <- c(length(annotation_EPIC$Relation_to_Island[which(annotation_EPIC$Relation_to_Island == "OpenSea")]), length(annotation_EPIC$Relation_to_Island[which(annotation_EPIC$Relation_to_Island == "Shore")]), length(annotation_EPIC$Relation_to_Island[which(annotation_EPIC$Relation_to_Island == "Shelf")]), length(annotation_EPIC$Relation_to_Island[which(annotation_EPIC$Relation_to_Island == "Island")]),
                      length(R_vs_No_HYPER$Relation_to_Island[which(R_vs_No_HYPER$Relation_to_Island == "OpenSea")]), length(R_vs_No_HYPER$Relation_to_Island[which(R_vs_No_HYPER$Relation_to_Island == "Shore")]), length(R_vs_No_HYPER$Relation_to_Island[which(R_vs_No_HYPER$Relation_to_Island == "Shelf")]), length(R_vs_No_HYPER$Relation_to_Island[which(R_vs_No_HYPER$Relation_to_Island == "Island")]),
                      length(R_vs_No_hypo$Relation_to_Island[which(R_vs_No_hypo$Relation_to_Island == "OpenSea")]), length(R_vs_No_hypo$Relation_to_Island[which(R_vs_No_hypo$Relation_to_Island == "Shore")]), length(R_vs_No_hypo$Relation_to_Island[which(R_vs_No_hypo$Relation_to_Island == "Shelf")]), length(R_vs_No_hypo$Relation_to_Island[which(R_vs_No_hypo$Relation_to_Island == "Island")]))

numbers_island_Recid <- c(length(annotation_EPIC$Relation_to_Island[which(annotation_EPIC$Relation_to_Island == "OpenSea")]), length(annotation_EPIC$Relation_to_Island[which(annotation_EPIC$Relation_to_Island == "Shore")]), length(annotation_EPIC$Relation_to_Island[which(annotation_EPIC$Relation_to_Island == "Shelf")]), length(annotation_EPIC$Relation_to_Island[which(annotation_EPIC$Relation_to_Island == "Island")]),
                          length(Recid_HYPER$Relation_to_Island[which(Recid_HYPER$Relation_to_Island == "OpenSea")]), length(Recid_HYPER$Relation_to_Island[which(Recid_HYPER$Relation_to_Island == "Shore")]), length(Recid_HYPER$Relation_to_Island[which(Recid_HYPER$Relation_to_Island == "Shelf")]), length(Recid_HYPER$Relation_to_Island[which(Recid_HYPER$Relation_to_Island == "Island")]),
                          length(Recid_hypo$Relation_to_Island[which(Recid_hypo$Relation_to_Island == "OpenSea")]), length(Recid_hypo$Relation_to_Island[which(Recid_hypo$Relation_to_Island == "Shore")]), length(Recid_hypo$Relation_to_Island[which(Recid_hypo$Relation_to_Island == "Shelf")]), length(Recid_hypo$Relation_to_Island[which(Recid_hypo$Relation_to_Island == "Island")]))

numbers_island_RECUR <- c(length(annotation_EPIC$Relation_to_Island[which(annotation_EPIC$Relation_to_Island == "OpenSea")]), length(annotation_EPIC$Relation_to_Island[which(annotation_EPIC$Relation_to_Island == "Shore")]), length(annotation_EPIC$Relation_to_Island[which(annotation_EPIC$Relation_to_Island == "Shelf")]), length(annotation_EPIC$Relation_to_Island[which(annotation_EPIC$Relation_to_Island == "Island")]),
                          length(RECUR_HYPER$Relation_to_Island[which(RECUR_HYPER$Relation_to_Island == "OpenSea")]), length(RECUR_HYPER$Relation_to_Island[which(RECUR_HYPER$Relation_to_Island == "Shore")]), length(RECUR_HYPER$Relation_to_Island[which(RECUR_HYPER$Relation_to_Island == "Shelf")]), length(RECUR_HYPER$Relation_to_Island[which(RECUR_HYPER$Relation_to_Island == "Island")]),
                          length(RECUR_hypo$Relation_to_Island[which(RECUR_hypo$Relation_to_Island == "OpenSea")]), length(RECUR_hypo$Relation_to_Island[which(RECUR_hypo$Relation_to_Island == "Shore")]), length(RECUR_hypo$Relation_to_Island[which(RECUR_hypo$Relation_to_Island == "Shelf")]), length(RECUR_hypo$Relation_to_Island[which(RECUR_hypo$Relation_to_Island == "Island")]))


numbers_island_RECUR_parcial <- c(length(annotation_EPIC$Relation_to_Island[which(annotation_EPIC$Relation_to_Island == "OpenSea")]), length(annotation_EPIC$Relation_to_Island[which(annotation_EPIC$Relation_to_Island == "Shore")]), length(annotation_EPIC$Relation_to_Island[which(annotation_EPIC$Relation_to_Island == "Shelf")]), length(annotation_EPIC$Relation_to_Island[which(annotation_EPIC$Relation_to_Island == "Island")]),
                                  length(RECUR_parcial_HYPER$Relation_to_Island[which(RECUR_parcial_HYPER$Relation_to_Island == "OpenSea")]), length(RECUR_parcial_HYPER$Relation_to_Island[which(RECUR_parcial_HYPER$Relation_to_Island == "Shore")]), length(RECUR_parcial_HYPER$Relation_to_Island[which(RECUR_parcial_HYPER$Relation_to_Island == "Shelf")]), length(RECUR_parcial_HYPER$Relation_to_Island[which(RECUR_parcial_HYPER$Relation_to_Island == "Island")]),
                                  length(RECUR_parcial_hypo$Relation_to_Island[which(RECUR_parcial_hypo$Relation_to_Island == "OpenSea")]), length(RECUR_parcial_hypo$Relation_to_Island[which(RECUR_parcial_hypo$Relation_to_Island == "Shore")]), length(RECUR_parcial_hypo$Relation_to_Island[which(RECUR_parcial_hypo$Relation_to_Island == "Shelf")]), length(RECUR_parcial_hypo$Relation_to_Island[which(RECUR_parcial_hypo$Relation_to_Island == "Island")]))

df_island_R <- data.frame(group_island_R,class_island,numbers_island_R)
df_island_R$class_island <- factor(df_island_R$class_island, levels = c("Island", "Shelf", "Shore", "OpenSea"))
df_island_R$group_island_R <- factor(df_island_R$group_island_R, levels = c("Back", "HYPER", "hypo"))

df_island_Recid <- data.frame(group_island_Recid,class_island,numbers_island_Recid)
df_island_Recid$class_island <- factor(df_island_Recid$class_island, levels = c("Island", "Shelf", "Shore", "OpenSea"))
df_island_Recid$group_island_Recid <- factor(df_island_Recid$group_island_Recid, levels = c("Back", "HYPER", "hypo"))

df_island_RECUR <- data.frame(group_island_RECUR,class_island,numbers_island_RECUR)
df_island_RECUR$class_island <- factor(df_island_RECUR$class_island, levels = c("Island", "Shelf", "Shore", "OpenSea"))
df_island_RECUR$group_island_RECUR <- factor(df_island_RECUR$group_island_RECUR, levels = c("Back", "HYPER", "hypo"))

df_island_RECUR_parcial <- data.frame(group_island_RECUR_parcial,class_island,numbers_island_RECUR_parcial)
df_island_RECUR_parcial$class_island <- factor(df_island_RECUR_parcial$class_island, levels = c("Island", "Shelf", "Shore", "OpenSea"))
df_island_RECUR_parcial$group_island_RECUR_parcial <- factor(df_island_RECUR_parcial$group_island_RECUR_parcial, levels = c("Back", "HYPER", "hypo"))


setwd(resultsDirectory)

pdf(file = "DMPs_genomic_distribution_CpGs_v2.pdf", height = 5, width = 4, onefile = T)
ggplot(df_island_R, aes(x = group_island_R)) + geom_bar(aes(weight= numbers_island_R, fill = class_island), position = 'fill') + scale_fill_manual(values = rev(c("#d09c9c", "#bf7272", "#a72b2b", "#761212"))) + theme_bw() + ggtitle(label = "R_vs_No \nCpG distribution CpGI status") + ylab(label = "Relative Frequency") + theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5)) + xlab("")
ggplot(df_island_Recid, aes(x = group_island_Recid)) + geom_bar(aes(weight= numbers_island_Recid, fill = class_island), position = 'fill') + scale_fill_manual(values = rev(c("#d09c9c", "#bf7272", "#a72b2b", "#761212"))) + theme_bw() + ggtitle(label = "Recid_vs_No \nCpG distribution CpGI status") + ylab(label = "Relative Frequency") + theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5)) + xlab("")
ggplot(df_island_RECUR, aes(x = group_island_RECUR)) + geom_bar(aes(weight= numbers_island_RECUR, fill = class_island), position = 'fill') + scale_fill_manual(values = rev(c("#d09c9c", "#bf7272", "#a72b2b", "#761212"))) + theme_bw() + ggtitle(label = "RECUR_vs_No \nCpG distribution CpGI status") + ylab(label = "Relative Frequency") + theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5)) + xlab("")
ggplot(df_island_RECUR_parcial, aes(x = group_island_RECUR_parcial)) + geom_bar(aes(weight= numbers_island_RECUR_parcial, fill = class_island), position = 'fill') + scale_fill_manual(values = rev(c("#d09c9c", "#bf7272", "#a72b2b", "#761212"))) + theme_bw() + ggtitle(label = "RECUR_vs_Parcial \nCpG distribution CpGI status") + ylab(label = "Relative Frequency") + theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5)) + xlab("")

#######################################  REGIONS  ############################################

group_region_R <- c(rep("Back",8), rep("HYPER", 8), rep("hypo", 8))
group_region_Recid <- c(rep("Back",8), rep("HYPER", 8), rep("hypo", 8))
group_region_RECUR <- c(rep("Back",8), rep("HYPER", 8), rep("hypo", 8))
group_region_RECUR_parcial <- c(rep("Back",8), rep("HYPER", 8), rep("hypo", 8))
class_region <- c(rep(c("Distal promoter", "Promoter (<=1kb)", "5' UTR", "Intron", "Exon", "3' UTR", "Downstream", "Distal Intergenic"), 3))

numbers_region_R <- c(length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Distal promoter")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Promoter (<=1kb)")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "5' UTR")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Intron")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Exon")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "3' UTR")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Downstream")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Distal Intergenic")]),
                      length(R_vs_No_HYPER$annotation[which(R_vs_No_HYPER$annotation == "Distal promoter")]), length(R_vs_No_HYPER$annotation[which(R_vs_No_HYPER$annotation == "Promoter (<=1kb)")]), length(R_vs_No_HYPER$annotation[which(R_vs_No_HYPER$annotation == "5' UTR")]), length(R_vs_No_HYPER$annotation[which(R_vs_No_HYPER$annotation == "Intron")]), length(R_vs_No_HYPER$annotation[which(R_vs_No_HYPER$annotation == "Exon")]), length(R_vs_No_HYPER$annotation[which(R_vs_No_HYPER$annotation == "3' UTR")]), length(R_vs_No_HYPER$annotation[which(R_vs_No_HYPER$annotation == "Downstream")]), length(R_vs_No_HYPER$annotation[which(R_vs_No_HYPER$annotation == "Distal Intergenic")]),
                      length(R_vs_No_hypo$annotation[which(R_vs_No_hypo$annotation == "Distal promoter")]), length(R_vs_No_hypo$annotation[which(R_vs_No_hypo$annotation == "Promoter (<=1kb)")]), length(R_vs_No_hypo$annotation[which(R_vs_No_hypo$annotation == "5' UTR")]), length(R_vs_No_hypo$annotation[which(R_vs_No_hypo$annotation == "Intron")]), length(R_vs_No_hypo$annotation[which(R_vs_No_hypo$annotation == "Exon")]), length(R_vs_No_hypo$annotation[which(R_vs_No_hypo$annotation == "3' UTR")]), length(R_vs_No_hypo$annotation[which(R_vs_No_hypo$annotation == "Downstream")]), length(R_vs_No_hypo$annotation[which(R_vs_No_hypo$annotation == "Distal Intergenic")]))

numbers_region_Recid <- c(length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Distal promoter")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Promoter (<=1kb)")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "5' UTR")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Intron")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Exon")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "3' UTR")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Downstream")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Distal Intergenic")]),
                          length(Recid_HYPER$annotation[which(Recid_HYPER$annotation == "Distal promoter")]), length(Recid_HYPER$annotation[which(Recid_HYPER$annotation == "Promoter (<=1kb)")]), length(Recid_HYPER$annotation[which(Recid_HYPER$annotation == "5' UTR")]), length(Recid_HYPER$annotation[which(Recid_HYPER$annotation == "Intron")]), length(Recid_HYPER$annotation[which(Recid_HYPER$annotation == "Exon")]), length(Recid_HYPER$annotation[which(Recid_HYPER$annotation == "3' UTR")]), length(Recid_HYPER$annotation[which(Recid_HYPER$annotation == "Downstream")]), length(Recid_HYPER$annotation[which(Recid_HYPER$annotation == "Distal Intergenic")]),
                          length(Recid_hypo$annotation[which(Recid_hypo$annotation == "Distal promoter")]), length(Recid_hypo$annotation[which(Recid_hypo$annotation == "Promoter (<=1kb)")]), length(Recid_hypo$annotation[which(Recid_hypo$annotation == "5' UTR")]), length(Recid_hypo$annotation[which(Recid_hypo$annotation == "Intron")]), length(Recid_hypo$annotation[which(Recid_hypo$annotation == "Exon")]), length(Recid_hypo$annotation[which(Recid_hypo$annotation == "3' UTR")]), length(Recid_hypo$annotation[which(Recid_hypo$annotation == "Downstream")]), length(Recid_hypo$annotation[which(Recid_hypo$annotation == "Distal Intergenic")]))

numbers_region_RECUR <- c(length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Distal promoter")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Promoter (<=1kb)")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "5' UTR")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Intron")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Exon")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "3' UTR")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Downstream")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Distal Intergenic")]),
                          length(RECUR_HYPER$annotation[which(RECUR_HYPER$annotation == "Distal promoter")]), length(RECUR_HYPER$annotation[which(RECUR_HYPER$annotation == "Promoter (<=1kb)")]), length(RECUR_HYPER$annotation[which(RECUR_HYPER$annotation == "5' UTR")]), length(RECUR_HYPER$annotation[which(RECUR_HYPER$annotation == "Intron")]), length(RECUR_HYPER$annotation[which(RECUR_HYPER$annotation == "Exon")]), length(RECUR_HYPER$annotation[which(RECUR_HYPER$annotation == "3' UTR")]), length(RECUR_HYPER$annotation[which(RECUR_HYPER$annotation == "Downstream")]), length(RECUR_HYPER$annotation[which(RECUR_HYPER$annotation == "Distal Intergenic")]),
                          length(RECUR_hypo$annotation[which(RECUR_hypo$annotation == "Distal promoter")]), length(RECUR_hypo$annotation[which(RECUR_hypo$annotation == "Promoter (<=1kb)")]), length(RECUR_hypo$annotation[which(RECUR_hypo$annotation == "5' UTR")]), length(RECUR_hypo$annotation[which(RECUR_hypo$annotation == "Intron")]), length(RECUR_hypo$annotation[which(RECUR_hypo$annotation == "Exon")]), length(RECUR_hypo$annotation[which(RECUR_hypo$annotation == "3' UTR")]), length(RECUR_hypo$annotation[which(RECUR_hypo$annotation == "Downstream")]), length(RECUR_hypo$annotation[which(RECUR_hypo$annotation == "Distal Intergenic")]))

numbers_region_RECUR_parcial <- c(length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Distal promoter")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Promoter (<=1kb)")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "5' UTR")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Intron")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Exon")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "3' UTR")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Downstream")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Distal Intergenic")]),
                                  length(RECUR_parcial_HYPER$annotation[which(RECUR_parcial_HYPER$annotation == "Distal promoter")]), length(RECUR_parcial_HYPER$annotation[which(RECUR_parcial_HYPER$annotation == "Promoter (<=1kb)")]), length(RECUR_parcial_HYPER$annotation[which(RECUR_parcial_HYPER$annotation == "5' UTR")]), length(RECUR_parcial_HYPER$annotation[which(RECUR_parcial_HYPER$annotation == "Intron")]), length(RECUR_parcial_HYPER$annotation[which(RECUR_parcial_HYPER$annotation == "Exon")]), length(RECUR_parcial_HYPER$annotation[which(RECUR_parcial_HYPER$annotation == "3' UTR")]), length(RECUR_parcial_HYPER$annotation[which(RECUR_parcial_HYPER$annotation == "Downstream")]), length(RECUR_parcial_HYPER$annotation[which(RECUR_parcial_HYPER$annotation == "Distal Intergenic")]),
                                  length(RECUR_parcial_hypo$annotation[which(RECUR_parcial_hypo$annotation == "Distal promoter")]), length(RECUR_parcial_hypo$annotation[which(RECUR_parcial_hypo$annotation == "Promoter (<=1kb)")]), length(RECUR_parcial_hypo$annotation[which(RECUR_parcial_hypo$annotation == "5' UTR")]), length(RECUR_parcial_hypo$annotation[which(RECUR_parcial_hypo$annotation == "Intron")]), length(RECUR_parcial_hypo$annotation[which(RECUR_parcial_hypo$annotation == "Exon")]), length(RECUR_parcial_hypo$annotation[which(RECUR_parcial_hypo$annotation == "3' UTR")]), length(RECUR_parcial_hypo$annotation[which(RECUR_parcial_hypo$annotation == "Downstream")]), length(RECUR_parcial_hypo$annotation[which(RECUR_parcial_hypo$annotation == "Distal Intergenic")]))


df_region_R <- data.frame(group_region_R,class_region,numbers_region_R)
df_region_R$class_region <- factor(df_region_R$class_region, levels = c("Distal promoter", "Promoter (<=1kb)", "5' UTR", "Intron", "Exon", "3' UTR", "Downstream", "Distal Intergenic"))
df_region_R$group_region_R <- factor(df_region_R$group_region_R, levels = c("Back", "HYPER", "hypo"))

df_region_Recid <- data.frame(group_region_Recid,class_region,numbers_region_Recid)
df_region_Recid$class_region <- factor(df_region_Recid$class_region, levels = c("Distal promoter", "Promoter (<=1kb)", "5' UTR", "Intron", "Exon", "3' UTR", "Downstream", "Distal Intergenic"))
df_region_Recid$group_region_Recid <- factor(df_region_Recid$group_region_Recid, levels = c("Back", "HYPER", "hypo"))

df_region_RECUR <- data.frame(group_region_RECUR,class_region,numbers_region_RECUR)
df_region_RECUR$class_region <- factor(df_region_RECUR$class_region, levels = c("Distal promoter", "Promoter (<=1kb)", "5' UTR", "Intron", "Exon", "3' UTR", "Downstream", "Distal Intergenic"))
df_region_RECUR$group_region_RECUR <- factor(df_region_RECUR$group_region_RECUR, levels = c("Back", "HYPER", "hypo"))

df_region_RECUR_parcial <- data.frame(group_region_RECUR_parcial,class_region,numbers_region_RECUR_parcial)
df_region_RECUR_parcial$class_region <- factor(df_region_RECUR_parcial$class_region, levels = c("Distal promoter", "Promoter (<=1kb)", "5' UTR", "Intron", "Exon", "3' UTR", "Downstream", "Distal Intergenic"))
df_region_RECUR_parcial$group_region_RECUR_parcial <- factor(df_region_RECUR_parcial$group_region_RECUR_parcial, levels = c("Back", "HYPER", "hypo"))


ggplot(df_region_R, aes(x = group_region_R)) + geom_bar(aes(weight= numbers_region_R, fill = class_region), position = 'fill') + scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")) + theme_bw() + ggtitle(label = "R_vs_No \nCpG distribution Region status") + ylab(label = "Relative Frequency") + theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5)) + xlab("")
ggplot(df_region_Recid, aes(x = group_region_Recid)) + geom_bar(aes(weight= numbers_region_Recid, fill = class_region), position = 'fill') + scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")) + theme_bw() + ggtitle(label = "Recid_vs_No \nCpG distribution Region status") + ylab(label = "Relative Frequency") + theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5)) + xlab("")
ggplot(df_region_RECUR, aes(x = group_region_RECUR)) + geom_bar(aes(weight= numbers_region_RECUR, fill = class_region), position = 'fill') + scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")) + theme_bw() + ggtitle(label = "RECUR_vs_No \nCpG distribution Region status") + ylab(label = "Relative Frequency") + theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5)) + xlab("")
ggplot(df_region_RECUR_parcial, aes(x = group_region_RECUR_parcial)) + geom_bar(aes(weight= numbers_region_RECUR_parcial, fill = class_region), position = 'fill') + scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")) + theme_bw() + ggtitle(label = "RECUR_vs_Parcial \nCpG distribution Region status") + ylab(label = "Relative Frequency") + theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5)) + xlab("")

dev.off()


######################################################################################################################
###########################################           OTHER GRAPHS          ##########################################
######################################################################################################################


# Establish the data frames and perform the plots
group_island_hyper <- c(rep("Back",4), rep("R_vs_No", 4), rep("Recid", 4), rep("RECUR", 4), rep("RECUR_parcial", 4))
group_island_hypo <- c(rep("Back",4), rep("R_vs_No", 4), rep("Recid", 4), rep("RECUR", 4), rep("RECUR_parcial", 4))
class_island_hyper <- c(rep(c("OpenSea", "Shore", "Shelf", "Island"), 5))
class_island_hypo <- c(rep(c("OpenSea", "Shore", "Shelf", "Island"), 5))
numbers_island_hyper <- c(length(annotation_EPIC$Relation_to_Island[which(annotation_EPIC$Relation_to_Island == "OpenSea")]), length(annotation_EPIC$Relation_to_Island[which(annotation_EPIC$Relation_to_Island == "Shore")]), length(annotation_EPIC$Relation_to_Island[which(annotation_EPIC$Relation_to_Island == "Shelf")]), length(annotation_EPIC$Relation_to_Island[which(annotation_EPIC$Relation_to_Island == "Island")]),
                          length(R_vs_No_HYPER$Relation_to_Island[which(R_vs_No_HYPER$Relation_to_Island == "OpenSea")]), length(R_vs_No_HYPER$Relation_to_Island[which(R_vs_No_HYPER$Relation_to_Island == "Shore")]), length(R_vs_No_HYPER$Relation_to_Island[which(R_vs_No_HYPER$Relation_to_Island == "Shelf")]), length(R_vs_No_HYPER$Relation_to_Island[which(R_vs_No_HYPER$Relation_to_Island == "Island")]),
                          length(Recid_HYPER$Relation_to_Island[which(Recid_HYPER$Relation_to_Island == "OpenSea")]), length(Recid_HYPER$Relation_to_Island[which(Recid_HYPER$Relation_to_Island == "Shore")]), length(Recid_HYPER$Relation_to_Island[which(Recid_HYPER$Relation_to_Island == "Shelf")]), length(Recid_HYPER$Relation_to_Island[which(Recid_HYPER$Relation_to_Island == "Island")]),
                          length(RECUR_HYPER$Relation_to_Island[which(RECUR_HYPER$Relation_to_Island == "OpenSea")]), length(RECUR_HYPER$Relation_to_Island[which(RECUR_HYPER$Relation_to_Island == "Shore")]), length(RECUR_HYPER$Relation_to_Island[which(RECUR_HYPER$Relation_to_Island == "Shelf")]), length(RECUR_HYPER$Relation_to_Island[which(RECUR_HYPER$Relation_to_Island == "Island")]),
                          length(RECUR_parcial_HYPER$Relation_to_Island[which(RECUR_parcial_HYPER$Relation_to_Island == "OpenSea")]), length(RECUR_parcial_HYPER$Relation_to_Island[which(RECUR_parcial_HYPER$Relation_to_Island == "Shore")]), length(RECUR_parcial_HYPER$Relation_to_Island[which(RECUR_parcial_HYPER$Relation_to_Island == "Shelf")]), length(RECUR_parcial_HYPER$Relation_to_Island[which(RECUR_parcial_HYPER$Relation_to_Island == "Island")]))

numbers_island_hypo <- c(length(annotation_EPIC$Relation_to_Island[which(annotation_EPIC$Relation_to_Island == "OpenSea")]), length(annotation_EPIC$Relation_to_Island[which(annotation_EPIC$Relation_to_Island == "Shore")]), length(annotation_EPIC$Relation_to_Island[which(annotation_EPIC$Relation_to_Island == "Shelf")]), length(annotation_EPIC$Relation_to_Island[which(annotation_EPIC$Relation_to_Island == "Island")]),
                         length(R_vs_No_Hypo$Relation_to_Island[which(R_vs_No_Hypo$Relation_to_Island == "OpenSea")]), length(R_vs_No_Hypo$Relation_to_Island[which(R_vs_No_Hypo$Relation_to_Island == "Shore")]), length(R_vs_No_Hypo$Relation_to_Island[which(R_vs_No_Hypo$Relation_to_Island == "Shelf")]), length(R_vs_No_Hypo$Relation_to_Island[which(R_vs_No_Hypo$Relation_to_Island == "Island")]),
                         length(Recid_Hypo$Relation_to_Island[which(Recid_Hypo$Relation_to_Island == "OpenSea")]), length(Recid_Hypo$Relation_to_Island[which(Recid_Hypo$Relation_to_Island == "Shore")]), length(Recid_Hypo$Relation_to_Island[which(Recid_Hypo$Relation_to_Island == "Shelf")]), length(Recid_Hypo$Relation_to_Island[which(Recid_Hypo$Relation_to_Island == "Island")]),
                         length(RECUR_Hypo$Relation_to_Island[which(RECUR_Hypo$Relation_to_Island == "OpenSea")]), length(RECUR_Hypo$Relation_to_Island[which(RECUR_Hypo$Relation_to_Island == "Shore")]), length(RECUR_Hypo$Relation_to_Island[which(RECUR_Hypo$Relation_to_Island == "Shelf")]), length(RECUR_Hypo$Relation_to_Island[which(RECUR_Hypo$Relation_to_Island == "Island")]),
                         length(RECUR_parcial_Hypo$Relation_to_Island[which(RECUR_parcial_Hypo$Relation_to_Island == "OpenSea")]), length(RECUR_parcial_Hypo$Relation_to_Island[which(RECUR_parcial_Hypo$Relation_to_Island == "Shore")]), length(RECUR_parcial_Hypo$Relation_to_Island[which(RECUR_parcial_Hypo$Relation_to_Island == "Shelf")]), length(RECUR_parcial_Hypo$Relation_to_Island[which(RECUR_parcial_Hypo$Relation_to_Island == "Island")]))

df_island_hyper <- data.frame(group_island_hyper,class_island_hyper,numbers_island_hyper)
df_island_hyper$class_island_hyper <- factor(df_island_hyper$class_island_hyper, levels = c("Island", "Shelf", "Shore", "OpenSea"))
df_island_hyper$group_island_hyper <- factor(df_island_hyper$group_island_hyper, levels = c("Back", "R_vs_No", "Recid", "RECUR", "RECUR_parcial"))

df_island_hypo <- data.frame(group_island_hypo,class_island_hypo,numbers_island_hypo)
df_island_hypo$class_island_hypo <- factor(df_island_hypo$class_island_hypo, levels = c("Island", "Shelf", "Shore", "OpenSea"))
df_island_hypo$group_island_hypo <- factor(df_island_hypo$group_island_hypo, levels = c("Back", "R_vs_No", "Recid", "RECUR", "RECUR_parcial"))

setwd(resultsDirectory)

pdf(file = "DMPs_genomic_distribution_CpGs.pdf", height = 5, width = 4, onefile = T)
ggplot(df_island_hyper, aes(x = group_island_hyper)) + geom_bar(aes(weight= numbers_island_hyper, fill = class_island_hyper), position = 'fill') + scale_fill_manual(values = rev(c("#d09c9c", "#bf7272", "#a72b2b", "#761212"))) + theme_bw() + ggtitle(label = "CpG Hyper distribution CpGI status") + ylab(label = "Relative Frequency") + theme(axis.text.x = element_text(angle = 90))
ggplot(df_island_hypo, aes(x = group_island_hypo)) + geom_bar(aes(weight= numbers_island_hypo, fill = class_island_hypo), position = 'fill') + scale_fill_manual(values = rev(c("#d09c9c", "#bf7272", "#a72b2b", "#761212"))) + theme_bw() + ggtitle(label = "CpG Hypo distribution CpGI status") + ylab(label = "Relative Frequency") + theme(axis.text.x = element_text(angle = 90))


###################################################################################

group_region_hyper <- c(rep("Back", 8), rep("R_vs_No", 8), rep("Recid", 8), rep("RECUR", 8), rep("RECUR_parcial", 8))
group_region_hypo <- c(rep("Back", 8), rep("R_vs_No", 8), rep("Recid", 8), rep("RECUR", 8), rep("RECUR_parcial", 8))
class_region_hyper <- c(rep(c("Distal promoter", "Promoter (<=1kb)", "5' UTR", "Intron", "Exon", "3' UTR", "Downstream", "Distal Intergenic"), 5))
class_region_hypo <- c(rep(c("Distal promoter", "Promoter (<=1kb)", "5' UTR", "Intron", "Exon", "3' UTR", "Downstream", "Distal Intergenic"), 5))

numbers_region_hyper <- c(length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Distal promoter")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Promoter (<=1kb)")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "5' UTR")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Intron")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Exon")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "3' UTR")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Downstream")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Distal Intergenic")]),
                          length(R_vs_No_HYPER$annotation[which(R_vs_No_HYPER$annotation == "Distal promoter")]), length(R_vs_No_HYPER$annotation[which(R_vs_No_HYPER$annotation == "Promoter (<=1kb)")]), length(R_vs_No_HYPER$annotation[which(R_vs_No_HYPER$annotation == "5' UTR")]), length(R_vs_No_HYPER$annotation[which(R_vs_No_HYPER$annotation == "Intron")]), length(R_vs_No_HYPER$annotation[which(R_vs_No_HYPER$annotation == "Exon")]), length(R_vs_No_HYPER$annotation[which(R_vs_No_HYPER$annotation == "3' UTR")]), length(R_vs_No_HYPER$annotation[which(R_vs_No_HYPER$annotation == "Downstream")]), length(R_vs_No_HYPER$annotation[which(R_vs_No_HYPER$annotation == "Distal Intergenic")]),
                          length(Recid_HYPER$annotation[which(Recid_HYPER$annotation == "Distal promoter")]), length(Recid_HYPER$annotation[which(Recid_HYPER$annotation == "Promoter (<=1kb)")]), length(Recid_HYPER$annotation[which(Recid_HYPER$annotation == "5' UTR")]), length(Recid_HYPER$annotation[which(Recid_HYPER$annotation == "Intron")]), length(Recid_HYPER$annotation[which(Recid_HYPER$annotation == "Exon")]), length(Recid_HYPER$annotation[which(Recid_HYPER$annotation == "3' UTR")]), length(Recid_HYPER$annotation[which(Recid_HYPER$annotation == "Downstream")]), length(Recid_HYPER$annotation[which(Recid_HYPER$annotation == "Distal Intergenic")]),
                          length(RECUR_HYPER$annotation[which(RECUR_HYPER$annotation == "Distal promoter")]), length(RECUR_HYPER$annotation[which(RECUR_HYPER$annotation == "Promoter (<=1kb)")]), length(RECUR_HYPER$annotation[which(RECUR_HYPER$annotation == "5' UTR")]), length(RECUR_HYPER$annotation[which(RECUR_HYPER$annotation == "Intron")]), length(RECUR_HYPER$annotation[which(RECUR_HYPER$annotation == "Exon")]), length(RECUR_HYPER$annotation[which(RECUR_HYPER$annotation == "3' UTR")]), length(RECUR_HYPER$annotation[which(RECUR_HYPER$annotation == "Downstream")]), length(RECUR_HYPER$annotation[which(RECUR_HYPER$annotation == "Distal Intergenic")]),
                          length(RECUR_parcial_HYPER$annotation[which(RECUR_parcial_HYPER$annotation == "Distal promoter")]), length(RECUR_parcial_HYPER$annotation[which(RECUR_parcial_HYPER$annotation == "Promoter (<=1kb)")]), length(RECUR_parcial_HYPER$annotation[which(RECUR_parcial_HYPER$annotation == "5' UTR")]), length(RECUR_parcial_HYPER$annotation[which(RECUR_parcial_HYPER$annotation == "Intron")]), length(RECUR_parcial_HYPER$annotation[which(RECUR_parcial_HYPER$annotation == "Exon")]), length(RECUR_parcial_HYPER$annotation[which(RECUR_parcial_HYPER$annotation == "3' UTR")]), length(RECUR_parcial_HYPER$annotation[which(RECUR_parcial_HYPER$annotation == "Downstream")]), length(RECUR_parcial_HYPER$annotation[which(RECUR_parcial_HYPER$annotation == "Distal Intergenic")]))

numbers_region_hypo <- c(length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Distal promoter")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Promoter (<=1kb)")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "5' UTR")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Intron")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Exon")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "3' UTR")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Downstream")]), length(annotation_EPIC$annotation[which(annotation_EPIC$annotation == "Distal Intergenic")]),
                         length(R_vs_No_Hypo$annotation[which(R_vs_No_Hypo$annotation == "Distal promoter")]), length(R_vs_No_Hypo$annotation[which(R_vs_No_Hypo$annotation == "Promoter (<=1kb)")]), length(R_vs_No_Hypo$annotation[which(R_vs_No_Hypo$annotation == "5' UTR")]), length(R_vs_No_Hypo$annotation[which(R_vs_No_Hypo$annotation == "Intron")]), length(R_vs_No_Hypo$annotation[which(R_vs_No_Hypo$annotation == "Exon")]), length(R_vs_No_Hypo$annotation[which(R_vs_No_Hypo$annotation == "3' UTR")]), length(R_vs_No_Hypo$annotation[which(R_vs_No_Hypo$annotation == "Downstream")]), length(R_vs_No_Hypo$annotation[which(R_vs_No_Hypo$annotation == "Distal Intergenic")]),
                         length(Recid_Hypo$annotation[which(Recid_Hypo$annotation == "Distal promoter")]), length(Recid_Hypo$annotation[which(Recid_Hypo$annotation == "Promoter (<=1kb)")]), length(Recid_Hypo$annotation[which(Recid_Hypo$annotation == "5' UTR")]), length(Recid_Hypo$annotation[which(Recid_Hypo$annotation == "Intron")]), length(Recid_Hypo$annotation[which(Recid_Hypo$annotation == "Exon")]), length(Recid_Hypo$annotation[which(Recid_Hypo$annotation == "3' UTR")]), length(Recid_Hypo$annotation[which(Recid_Hypo$annotation == "Downstream")]), length(Recid_Hypo$annotation[which(Recid_Hypo$annotation == "Distal Intergenic")]),
                         length(RECUR_Hypo$annotation[which(RECUR_Hypo$annotation == "Distal promoter")]), length(RECUR_Hypo$annotation[which(RECUR_Hypo$annotation == "Promoter (<=1kb)")]), length(RECUR_Hypo$annotation[which(RECUR_Hypo$annotation == "5' UTR")]), length(RECUR_Hypo$annotation[which(RECUR_Hypo$annotation == "Intron")]), length(RECUR_Hypo$annotation[which(RECUR_Hypo$annotation == "Exon")]), length(RECUR_Hypo$annotation[which(RECUR_Hypo$annotation == "3' UTR")]), length(RECUR_Hypo$annotation[which(RECUR_Hypo$annotation == "Downstream")]), length(RECUR_Hypo$annotation[which(RECUR_Hypo$annotation == "Distal Intergenic")]),
                         length(RECUR_parcial_Hypo$annotation[which(RECUR_parcial_Hypo$annotation == "Distal promoter")]), length(RECUR_parcial_Hypo$annotation[which(RECUR_parcial_Hypo$annotation == "Promoter (<=1kb)")]), length(RECUR_parcial_Hypo$annotation[which(RECUR_parcial_Hypo$annotation == "5' UTR")]), length(RECUR_parcial_Hypo$annotation[which(RECUR_parcial_Hypo$annotation == "Intron")]), length(RECUR_parcial_Hypo$annotation[which(RECUR_parcial_Hypo$annotation == "Exon")]), length(RECUR_parcial_Hypo$annotation[which(RECUR_parcial_Hypo$annotation == "3' UTR")]), length(RECUR_parcial_Hypo$annotation[which(RECUR_parcial_Hypo$annotation == "Downstream")]), length(RECUR_parcial_Hypo$annotation[which(RECUR_parcial_Hypo$annotation == "Distal Intergenic")]))

                      
df_region_hyper <- data.frame(group_region_hyper,class_region_hyper,numbers_region_hyper)
df_region_hyper$class_region_hyper <- factor(df_region_hyper$class_region_hyper, levels = c("Distal promoter", "Promoter (<=1kb)", "5' UTR", "Intron", "Exon", "3' UTR", "Downstream", "Distal Intergenic")) 
df_region_hyper$group_region_hyper <- factor(df_region_hyper$group_region_hyper, levels = c("Back", "R_vs_No", "Recid", "RECUR", "RECUR_parcial"))

df_region_hypo <- data.frame(group_region_hypo,class_region_hypo,numbers_region_hypo)
df_region_hypo$class_region_hypo <- factor(df_region_hypo$class_region_hypo, levels = c("Distal promoter", "Promoter (<=1kb)", "5' UTR", "Intron", "Exon", "3' UTR", "Downstream", "Distal Intergenic"))  
df_region_hypo$group_region_hypo <- factor(df_region_hypo$group_region_hypo, levels = c("Back", "R_vs_No", "Recid", "RECUR", "RECUR_parcial"))

ggplot(df_region_hyper, aes(x = group_region_hyper)) + geom_bar(aes(weight = numbers_region_hyper, fill = class_region_hyper), position = 'fill') + scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")) + theme_bw() + ggtitle(label = "CpG Hyper distribution Region status") + ylab(label = "Relative Frequency") + theme(axis.text.x = element_text(angle = 90))
ggplot(df_region_hypo, aes(x = group_region_hypo)) + geom_bar(aes(weight = numbers_region_hypo, fill = class_region_hypo), position = 'fill') + scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")) + theme_bw() + ggtitle(label = "CpG Hypo distribution Region status") + ylab(label = "Relative Frequency") + theme(axis.text.x = element_text(angle = 90))

dev.off()
