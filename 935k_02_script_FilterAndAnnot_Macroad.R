#######################################################################
# PROJECT: Non-functioning Pituitary Macroadenomas
#######################################################################
### Script for the analysis of Methylation array V2 data (EPICv2 935K)
### 
### Script for the loading, filtering and semi-exploratory analysis of 
### Methylation array data (935k v2) (based on SeSaMe)        
###
### Author: Rocío G. Urdinguio
########################################################################

###################################################################
# Libraries
###################################################################

# These are V2 EPIC arrays! We will use jokergoo's custom annotations (thanks, bioinformagician)

#BiocManager::install("jokergoo/IlluminaHumanMethylationEPICv2manifest")
#BiocManager::install("jokergoo/IlluminaHumanMethylationEPICv2anno.20a1.hg38")
# install.packages('R.utils')

library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(ggfortify)
library(stringr)
library(gridExtra)
library(grid)
library(readxl)
library(patchwork)
library(qs)
library(data.table)

library(sesame)
library(sesameData)

sesameDataCache()

###################################################################
# Sources 
###################################################################

library(tidyverse)
getCurrentFileLocation <- function() # Adapted from https://stackoverflow.com/questions/47044068/get-the-path-of-current-script
{
  this_file <- commandArgs() %>% 
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file)==0)
  {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  
  return(this_file)
  
}

source(file.path(dirname(getCurrentFileLocation()),"05_Script_Annotation_after_Sesame_processing.R"))

script_prefix <- stringr::str_split_fixed(basename(getCurrentFileLocation()),"_",2)[,1]
script_name <- tools::file_path_sans_ext(basename(getCurrentFileLocation()))

dir_output_script <- file.path(dir_output,script_name)
dir.create(dir_output_script, recursive = T)

#################################################################
### SET UP THE PATHS AND DIRECTORIES
#################################################################

# Set up basedir
getwd()
setwd("..")
getwd()
basedir <- getwd()

#IDATs localization:
idat_dir <- file.path(basedir, "Raw_Data_935K/")

# Set up the directory for the results:
results_dir <- file.path(basedir, "PDF/")
tables_dir <- file.path(basedir, "Tables/")
SS_dir <- file.path(basedir, "SampleSheet/")
databases_dir <- file.path(basedir, "Databases/")

##################################################################
### LOAD BETAS AND phenodata
##################################################################

# Saved as "betas_prep_collapsed.txt" at the script: 01_Rmd_Sesame_Macroad935_v2.Rmd
#betas_prep_collapsed <- openSesame(idat_dir, BPPARAM = BiocParallel::SnowParam(4), 
#                                  func = getBetas, # specify func = NULL to get a SigDF list of objects
#                                  prep="QCDPB", # recommended for EPICv2
#                                  collapseToPfx = TRUE,  #for EPICv2 probes: collapse names with suffixes. The default method is mean between replicate probes
#                                  collapseMethod = "mean",
#                                  mask = TRUE) # Raul: enable masking out probes)

betas_prep_collapsed <- read.delim(file = paste(tables_dir, "betas_prep_collapsed.txt", 
                                                sep = "/"))
betas_prep_collapsed <- as.matrix(betas_prep_collapsed)

phenodata <- fread(file.path(tables_dir, "phenodata.txt"))

#Clean and prepare phenodata table:
phenodata <- read.csv(file.path(SS_dir, "Samplesheet_Macroadenomas_935k_phenodata.csv"),
                      sep = ",")
phenodata[,"Sample_Name"] <- NULL
phenodata[,"Sample_Well"] <- NULL
phenodata[,"Sample_Plate"] <- NULL
phenodata[,"GS_ID"] <- NULL
phenodata[,"CustomerID"] <- NULL
phenodata[,"Orig_Sample_Code"] <- NULL
phenodata[,"Other_IDs"] <- NULL
phenodata[,"ID_Cris_Alicante"] <- NULL
phenodata[,"ID_Cris_Sevilla"] <- NULL
phenodata[,"Sample_Type"] <- NULL

phenodata[,"Birth_Date"] <- NULL
phenodata[,"Diagnosis_Date"] <- NULL
phenodata[,"Tumoral_type"] <- NULL
phenodata[,"Recurrence"] <- NULL
phenodata[,"Relapse"] <- NULL
phenodata[,"Sample_Group"] <- NULL
phenodata[,"Recur_Dates"] <- NULL
phenodata[,"Recur1_date"] <- NULL
phenodata[,"Recur2_date"] <- NULL
phenodata[,"Relap_Dates"] <- NULL
phenodata[,"Relap1_date"] <- NULL
phenodata[,"Relap2_date"] <- NULL

phenodata[,"Ki67_percentage"] <- NULL
phenodata[,"Comorbidities_all"] <- NULL
phenodata[,"FSH"] <- NULL
phenodata[,"LH"] <- NULL
phenodata[,"GH"] <- NULL
phenodata[,"GH_range"] <- NULL
phenodata[,"ACTH"] <- NULL
phenodata[,"ACTH_range"] <- NULL
phenodata[,"PRL"] <- NULL
phenodata[,"PRL_range"] <- NULL
phenodata[,"TSH"] <- NULL
phenodata[,"PMR"] <- NULL

#As factor
phenodata$Sentrix_ID <- as.factor(phenodata$Sentrix_ID)
phenodata$Hospital <- as.factor(phenodata$Hospital)
phenodata$Sex <- as.factor(phenodata$Sex)
phenodata$Diagnosis_Age_Interval <- as.factor(phenodata$Diagnosis_Age_Interval)
phenodata$Histological_Subclassification <- as.factor(phenodata$Histological_Subclassification)
phenodata$Recur_Relap <- as.factor(phenodata$Recur_Relap)
phenodata$Group <- as.factor(phenodata$Group)
phenodata$Invasion_grade_Knosp_scale <- as.factor(phenodata$Invasion_grade_Knosp_scale)
phenodata$Cavernous_sinus_invasion <- as.factor(phenodata$Cavernous_sinus_invasion)
phenodata$Resection_grade <- as.factor(phenodata$Resection_grade)
phenodata$Comorb_HTA <- as.factor(phenodata$Comorb_HTA)
phenodata$Comorb_Hypopituitarism <- as.factor(phenodata$Comorb_Hypopituitarism)
phenodata$Comorb_Hypothyroidism <- as.factor(phenodata$Comorb_Hypothyroidism)
phenodata$Comorb_Diabetes <- as.factor(phenodata$Comorb_Diabetes)
phenodata$Comorb_other <- as.factor(phenodata$Comorb_other)
phenodata$Smoking <- as.factor(phenodata$Smoking)
phenodata$Deceased <- as.factor(phenodata$Deceased) # Only one level
phenodata$Transplant_or_Ttm <- as.factor(phenodata$Transplant_or_Ttm) # Only one level
phenodata$FSH_range <- as.factor(phenodata$FSH_range)
phenodata$LH_range <- as.factor(phenodata$LH_range)
phenodata$TSH_range <- as.factor(phenodata$TSH_range)
phenodata$IHQ_FT <- as.factor(phenodata$IHQ_FT) # Only one level
phenodata$Array_version <- as.factor(phenodata$Array_version) # Only one level

#Asign names to the columns in bvals
phenodata$Sentrix_ID <- paste("X", phenodata$Sentrix_ID, sep = "")
phenodata$Sentrix_IDPos <- paste(phenodata$Sentrix_ID, phenodata$Sentrix_Position, sep = "_")
names <- data.frame(phenodata$Sentrix_IDPos, phenodata$Customer_M_ID)

betas <- as.data.frame(betas_prep_collapsed)
betas <- betas_prep_collapsed[complete.cases(betas), ]

names_ordered <- names[order(match(phenodata[,"Sentrix_IDPos"], colnames(betas))),]
colnames(betas) <- names_ordered$phenodata.Customer_M_ID


###################################################################
# Load annotations
###################################################################

########################################
# Load EPICv2 annotations from sesame
########################################

# From: https://doi.org/10.1186/s43682-023-00021-5, https://zwdzwd.github.io/InfiniumAnnotation, https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPICv2/

# Manifest
annotationEPIC <- fread(file.path(databases_dir,"EPICv2.hg38.manifest.tsv.gz"))

# Genes
annotationEPIC.genes <- fread(file.path(databases_dir,"EPICv2.hg38.manifest.gencode.v41.tsv.gz"))

# SNP, alignment mask
annotationEPIC.mask <- fread(file.path(databases_dir,"EPICv2.hg38.mask.tsv.gz"))

# SNP information
annotationEPIC.snps <- fread(file.path(databases_dir,"EPICv2.hg38.snp.tsv.gz"))

# Other (CGIs, etc.)
annotationEPIC.cgi <- fread(file.path(databases_dir,"CGI.20220904.gz"))

########################################
# Remove replicate probes to have same probes in annotation and data
########################################

annotationEPIC$ID2 <- str_split_fixed(annotationEPIC$Probe_ID,"_",2)[,1]

annotationEPIC <- annotationEPIC[!duplicated(annotationEPIC$ID2),]

table(rownames(betas) == annotationEPIC$ID2)

###################################################################
# Filter probes
###################################################################

# The betas come masked by sesame because we used mask = TRUE
table(!complete.cases(betas))

# Not all masked probes are masked in all samples
table(apply(betas,1,function(x){sum(is.na(x)) == length(x)}))
table(apply(betas,1,function(x){sum(is.na(x)) > 5}))

toremove <- list()

########################################
# sesame pipeline mask: detection p-values, mapping probes, SNPs
########################################

# Removing probes with NA in at least 1 sample
toremove$sesame.mask <- rownames(betas)[!complete.cases(betas)]

########################################
# non-CpG probes
########################################

# The new EPICv2 arrays have different types of probes
# cg: C from CpG
# ch: C from non-CpG
# ct: control
# nv: mutations
# rs: SNP
table(substr(rownames(betas),1,2))

noncg <- c("ch","ct","nt","rs","nv")
toremove$toremove.noncg <- rownames(betas)[substr(rownames(betas),1,2) %in% noncg]

########################################
# Replicate probes (EPIC V2)
########################################

# The V2 also includes replicate probes. Most replicates have a ~1 correlation (see Fig. S2D-F in https://epicom.biomedcentral.com/articles/10.1186/s43682-023-00021-5)

# Removal already done by sesame because we used collapseToPfx = TRUE: collapse names with suffixes. The default method is mean between replicate probes

########################################
# Sex chromosome probes
########################################

toremove$toremove.chrsex <- annotationEPIC$ID2[grepl("chrX|chrY",annotationEPIC$CpG_chrm)]
  
########################################
# Other chromosome probes
########################################

# Remove some probes with NA chr and other chrs
toremove$toremove.chrother <- annotationEPIC$ID2[grepl("_alt|chrM",annotationEPIC$CpG_chrm) | is.na(annotationEPIC$CpG_chrm)]

########################################
# Cross-reactive, multimapping probes: EPIC v2
########################################

# From this paper: https://doi.org/10.1186/s12864-024-10027-5
cross.peters <- fread(file.path(databases_dir,"12864_2024_10027_MOESM8_ESM.csv"))

toremove$toremove.crosspeters <- str_split_fixed(unique(cross.peters$ProbeID),"_",2)[,1]

########################################
# Cross-reactive, multimapping probes: EPIC v1
########################################

# The new EPIC V2 version has removed and cleaned a lot of cross-reactive, multimapping probes. But there are some left.

# Multimapping/non-specific probes discovered in 450K array (2013, Chen et al., Epigenetics)
# File in: https://github.com/Jfortin1/funnorm_repro/blob/master/bad_probes/48639-non-specific-probes-Illumina450k.xlsx

multimapChen450K.cpg <- read_excel(path=file.path(databases_dir,"48639-non-specific-probes-Illumina450k.xlsx"), sheet = 1)
multimapChen450K.chh <- read_excel(path=file.path(databases_dir,"48639-non-specific-probes-Illumina450k.xlsx"), sheet = 2)

# Cross-reactive probes discovered in EPIC array (2016, Pidsley et al.)
# File in: https://github.com/markgene/maxprobes/blob/master/inst/extdata/13059_2016_1066_MOESM1_ESM.csv

cross.reactiveEPIC <- read.csv(file=file.path(databases_dir,"13059_2016_1066_MOESM1_ESM.csv"),header=T,sep=",")

# To remove
toremove$toremove.epicv1 <- unique(c(as.character(multimapChen450K.cpg$TargetID),as.character(multimapChen450K.chh$TargetID),as.character(cross.reactiveEPIC$X)))

# These other are redundant and not used:

# Cross-reactive probes discovered in EPIC array for CpG and CH probes (2016, McCartney DL et al.)

# setwd("/media/rauldiul/Externo/2021_PROYECTO_FCIEN/other/meth_filters")
# 
# cross.reactiveEPIC2 <- read.table(file="mmc2.txt",sep="\t",header=F) # These are mostly the same as cross.reactiveEPIC
# cross.reactiveEPIC2.ch <- read.table(file="mmc3.txt",sep="\t",header=F)
# 
# cross.reactiveEPIC2.probes <- cross.reactiveEPIC2$V1 # They are all in cross.reactiveEPIC.probes
# cross.reactiveEPIC2.ch.probes <- cross.reactiveEPIC2.ch$V1 # They are all in multimapChen450K.chh.probes
# length(intersect(cross.reactiveEPIC2.probes,cross.reactiveEPIC.probes))
# length(intersect(cross.reactiveEPIC2.ch.probes,multimapChen450K.chh.probes))

########################################
# sesame general mask
########################################

# This is the recommended mask from the sesame EPICv2 paper: https://doi.org/10.1186/s43682-023-00021-5, also described here https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10916044/
# See header description in https://zwdzwd.github.io/InfiniumAnnotation/mask2.html: the M_General column indicates the recommended masking
toremove$toremove.sesame.generalmask <- str_split_fixed(annotationEPIC.mask$Probe_ID[annotationEPIC.mask$M_general == TRUE],"_",2)[,1]

########################################
# Other sesame EPICv2 filters: SNP, mapping
########################################

# Sesame internally has some masks related to probe quality mapping and SNPs:
listAvailableMasks("EPICv2")
sesameDataList()

# We get them like this:
mask.general <- sesameDataGet("KYCG.EPICv2.Mask.20230314")

# These are actually already filtered because we use the prep="QCDPB" pipeline, which includes "Q" filtering

# But if we want to explore them:

  # These are contained in the external mask information obtained from https://zwdzwd.github.io/InfiniumAnnotation, the recommended mask from the sesame EPICv2 paper: https://doi.org/10.1186/s43682-023-00021-5, also described here https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10916044/
  length(unique(unlist(mask.general)))
  table(unique(unlist(mask.general)) %in% annotationEPIC.mask$Probe_ID)
  
  # The external mask information recommends filtering certain probes (See header description in https://zwdzwd.github.io/InfiniumAnnotation/mask2.html: the M_General column indicates the recommended masking)
  mask.general.external <- unique(annotationEPIC.mask$Probe_ID[annotationEPIC.mask$M_general == TRUE])
  
  # This should be equivalent to doing this:
  mask.general2 <- unique(unlist(mask.general[c("M_mapping","M_nonuniq","M_SNPcommon_5pt","M_1baseSwitchSNPcommon_5pt","M_2extBase_SNP_common_5pt")]))
  
  # It is similar but not the same
  # We will use the external recommended mask, because it is a bit bigger and contains the other
  length(mask.general.external); length(mask.general2)
  length(intersect(mask.general.external,mask.general2))
  
  # We can see that these probes are already masked with NA in our data because we used the prep="QCDPB" pipeline, which includes "Q" filtering
  foo <- betas[complete.cases(betas),]
  
  table(str_split_fixed(mask.general2,"_",2)[,1] %in% rownames(foo))
  table(str_split_fixed(mask.general.external,"_",2)[,1] %in% rownames(foo))

# So we include them to be removed but they are already removed

toremove$toremove.sesame.generalmask <- str_split_fixed(annotationEPIC.mask$Probe_ID[annotationEPIC.mask$M_general == TRUE],"_",2)[,1]
  
########################################
# sesame SNP information
########################################

# This is not a mask I think, it is just information on the SNPs
# View(annotationEPIC.snps)

########################################
# Other multimapping probes
########################################

# I'M CURRENTLY NOT USING THESE BECAUSE THEY ARE NOT PUBLISHED ANYWHERE

# # Exclude multi-mapping probes:  450K, Obtained from https://github.com/sirselim/illumina450k_filtering
# multimapProbes <- read.csv(file=paste(dataDirectory, "HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt", sep="/"), header=FALSE, as.is=TRUE)

########################################
# Final probes to remove
########################################

toremove.final <- unique(unlist(toremove))

betas.f <- betas[!(rownames(betas) %in% toremove.final |
                             str_split_fixed(rownames(betas),"_",2)[,1] %in% toremove.final),]

dim(betas);dim(betas.f)

###################################################################
# Expand sesame's array annotation: compare to Peters manifest
###################################################################

########################################
# Recover strand information
########################################

# Using the gene annotation to get the strand
annotationEPIC$probe_strand <- annotationEPIC.genes$probe_strand[match(annotationEPIC$Probe_ID,annotationEPIC.genes$probeID)]

# Note that this strand information is not the same as the "T" and "B" suffix from the probe names
table(substr(str_split_fixed(annotationEPIC$Probe_ID,"_",2)[,2],1,1), annotationEPIC$probe_strand)

# We do a sanity check with the prior EPICv1 annotation and find that the strand annotations are OPPOSITE between the 2 arrays!
#BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annotationEPICv1 <- as.data.frame(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))

# Between Sesame and Peters they do agree...
inters <- intersect(annotationEPICv1$Name,annotationEPIC$ID2)
table(annotationEPICv1[inters,]$strand == annotationEPIC$probe_strand[match(inters,annotationEPIC$ID2)])

########################################
# Exploring the improved manifest from Peters et al: 10.1186/s12864-024-10027-5
########################################

# I don't particularly care for sesame's array annotation because it is missing things such as strand, and the CpG coordinates are tricky (>1 bp)
# Thus, we will use the expanded manifest from Peters et al

# Loading the improved manifest from Peters et al: 10.1186/s12864-024-10027-5
# From: https://bioconductor.org/packages/release/data/annotation/manuals/EPICv2manifest/man/EPICv2manifest.pdf
annotationEPIC.peters <- fread(file.path(databases_dir,"12864_2024_10027_MOESM4_ESM.csv"))

# Build new strand variable
library(dplyr)
annotationEPIC.peters$strand <- case_when(
  annotationEPIC.peters$Strand_FR == "F" ~ "+",
  annotationEPIC.peters$Strand_FR == "R" ~ "-",
  annotationEPIC.peters$Strand_FR == "0" ~ NA
  )

# Similar results to sesame's strand annotation
table(annotationEPIC.peters$strand == annotationEPIC$probe_strand[match(annotationEPIC.peters$IlmnID,annotationEPIC$Probe_ID)])

# Improve annotations: select single priority for CpG islands
annotationEPIC.peters$Relation_to_UCSC_CpG_Island2 <- sapply(annotationEPIC.peters$Relation_to_UCSC_CpG_Island,function(x){
  
  final <- "OpenSea"
  
  for(location in c("S_Shelf","N_Shelf","S_Shore","N_Shore","Island")){ # using this priority
    if(grepl(location,x)){
      final <- location
    }
  }
  
  final
  
})

# Similar results to sesame's island annotation
table(gsub("CGI;","",annotationEPIC.cgi$Knowledgebase) == gsub("N_|S_","",annotationEPIC.peters$Relation_to_UCSC_CpG_Island2[match(annotationEPIC.cgi$Probe_ID,annotationEPIC.peters$IlmnID)]))

# Also, the flagged cross-hybridization probes which we filtered out previously (see toremove$toremove.crosspeters)are also included in the other peters et al. data.frame
length(annotationEPIC.peters$IlmnID[annotationEPIC.peters$CH_WGBS_evidence == "Y"])
table(annotationEPIC.peters$IlmnID[annotationEPIC.peters$CH_WGBS_evidence == "Y"] %in% unique(cross.peters$ProbeID))

# With regards to gene annotations, I also do not like this manifest because the gene annotations are only these categories:
# so we know nothing about intronic probes, and the paper is not too clear on it
sort(table(str_split_fixed(annotationEPIC.peters$UCSC_RefGene_Group,";",2)[,1]))

###################################################################
# Expand sesame's array annotation: annotate to CpG islands and genes
###################################################################

########################################
# For CpG islands we can use sesame's annotation
########################################

annotationEPIC$Relation_to_CpG_Island <- annotationEPIC.cgi$Knowledgebase[match(annotationEPIC$Probe_ID,annotationEPIC.cgi$Probe_ID)]

########################################
# For genes we will do a custom annotation
########################################

#BiocManager::install("ChIPseeker")
library(ChIPseeker)
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

annotationEPIC2 <- annotationEPIC[!is.na(annotationEPIC$CpG_chrm) & !is.na(annotationEPIC$CpG_beg) & !is.na(annotationEPIC$CpG_end),c("CpG_chrm","CpG_beg","CpG_end","probe_strand","Probe_ID")]

annotationEPIC.gr <- makeGRangesFromDataFrame(annotationEPIC2,
                                              keep.extra.columns = T,
                                              seqnames.field = "CpG_chrm",
                                              start.field = "CpG_beg",
                                              end.field = "CpG_end",
                                              strand.field = "probe_strand",
                                              starts.in.df.are.0based = TRUE)

anno.chip <- annotatePeak(annotationEPIC.gr,
                     tssRegion = c(-1000, 0),
                     TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                     level = "transcript", # Annotate at the transcript level: we care for TSS of transcripts
                     assignGenomicAnnotation = TRUE,
                     genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                                                   "Downstream", "Intergenic"),
                     annoDb = NULL,
                     addFlankGeneInfo = FALSE,
                     flankDistance = 5000,
                     sameStrand = FALSE, # take into account the strand of the gene annotation
                     ignoreOverlap = FALSE,
                     ignoreUpstream = FALSE,
                     ignoreDownstream = FALSE,
                     overlap = "TSS",
                     verbose = TRUE
)

anno.chip <- as.data.frame(anno.chip)
colnames(anno.chip) <- paste0("chipseeker_",colnames(anno.chip))

annotationEPIC <- cbind(annotationEPIC,anno.chip[match(annotationEPIC$Probe_ID,anno.chip$chipseeker_Probe_ID),7:15])

# Simplify annotation
annotationEPIC$chipseeker_annotation2 <- case_when(
  grepl("Exon",annotationEPIC$chipseeker_annotation, ignore.case = T) ~ "Exon",
  grepl("Intron",annotationEPIC$chipseeker_annotation, ignore.case = T) ~ "Intron",
  grepl("Intergenic",annotationEPIC$chipseeker_annotation, ignore.case = T) ~ "Intergenic",
  grepl("Downstream",annotationEPIC$chipseeker_annotation, ignore.case = T) ~ "Downstream",
  grepl("Promoter",annotationEPIC$chipseeker_annotation, ignore.case = T) ~ "Promoter",
  grepl("3' UTR",annotationEPIC$chipseeker_annotation, ignore.case = T) ~ "3UTR",
  grepl("5' UTR",annotationEPIC$chipseeker_annotation, ignore.case = T) ~ "5UTR",
  T ~ annotationEPIC$chipseeker_annotation
)

###################################################################
# Data exploration
###################################################################

########################################
# Infer sex
########################################

#BiocManager::install("wateRmelon")
library(wateRmelon)
sex <- wateRmelon::estimateSex(betas)

phenodata$sexPred_predictedSex <- ifelse(sex$predicted_sex[match(phenodata$Customer_M_ID,rownames(sex))] == "Male","H","M")

########################################
# QC explore
########################################

pdf(file = file.path(results_dir, paste0(script_prefix,"_","QC_explore.pdf")), height = 6.5, width = 6.5)

layout(matrix(1:3,nrow=3,ncol=1))
# Check histograms and boxplots

hist(as.matrix(betas),breaks=1000, main = "noob")
hist(as.matrix(betas.f),breaks=1000, main = "noob, filtered") 

boxplot(as.matrix(betas), main = "noob")
boxplot(as.matrix(betas.f), main = "noob, filtered")

densityPlot(as.matrix(betas), main = "noob", legend=FALSE)
densityPlot(as.matrix(betas.f), main = "noob, filtered", legend=FALSE)

# Check data (prcomp will fail with NAs, use raw.f)
pca.raw <- prcomp(t(betas[complete.cases(betas),]))
pca.f <- prcomp(t(betas.f))

tidy.raw <- data.frame(
  Dim1 = pca.raw$x[,1],
  Dim2 = pca.raw$x[,2],
  ID = phenodata$Customer_M_ID,
  sex = phenodata$sexPred_predictedSex)

tidy.f <- data.frame(
  Dim1 = pca.f$x[,1],
  Dim2 = pca.f$x[,2],
  ID = phenodata$Customer_M_ID,
  sex = phenodata$sexPred_predictedSex)

varexp.raw <- (pca.raw$sdev^2 / sum(pca.raw$sdev^2)) * 100
varexp.f <- (pca.f$sdev^2 / sum(pca.f$sdev^2)) * 100

dimX <- 1
dimY <- 2

#install.packages('ggthemes')
#(ggthemes)
# For each variable
p1 <- ggplot(data = data.frame(x = pca.raw$x[,dimX], y = pca.raw$x[,dimY]), aes(y = y, x = x, label = phenodata$Customer_M_ID)) +
  geom_point(shape = 19, size = 1) + geom_text(hjust=0, vjust=0) +
  xlab(paste0("Dim",dimX," (var: ",round(varexp.raw[dimX],1),"%)")) + ylab(paste0("Dim",dimY," (var: ",round(varexp.raw[dimY],1),"%)")) +
  ggtitle(paste0("PCA raw")) +
  theme_minimal() #ggtheme 

p2 <- ggplot(data = data.frame(x = pca.f$x[,dimX], y = pca.f$x[,dimY]), aes(y = y, x = x, label = phenodata$Customer_M_ID)) +
  geom_point(shape = 19, size = 1) + geom_text(hjust=0, vjust=0) +
  xlab(paste0("Dim",dimX," (var: ",round(varexp.raw[dimX],1),"%)")) + ylab(paste0("Dim",dimY," (var: ",round(varexp.raw[dimY],1),"%)")) +
  ggtitle(paste0("PCA noob")) +
  theme_minimal() #ggtheme 

p3 <- ggplot(data = data.frame(x = pca.raw$x[,dimX], y = pca.raw$x[,dimY]), aes(y = y, x = x, label = phenodata$Customer_M_ID, color = phenodata$sexPred_predictedSex)) +
  geom_point(shape = 19, size = 1) + geom_text(vjust="inward",hjust="inward") +
  xlab(paste0("Dim",dimX," (var: ",round(varexp.raw[dimX],1),"%)")) + ylab(paste0("Dim",dimY," (var: ",round(varexp.raw[dimY],1),"%)")) +
  ggtitle(paste0("PCA raw")) +
  scale_color_discrete(name = "Sex") +
  theme_minimal() #ggtheme 

p4 <- ggplot(data = data.frame(x = pca.f$x[,dimX], y = pca.f$x[,dimY]), aes(y = y, x = x, label = phenodata$Customer_M_ID, color = phenodata$sexPred_predictedSex)) +
  geom_point(shape = 19, size = 1) + geom_text(vjust="inward",hjust="inward") +
  xlab(paste0("Dim",dimX," (var: ",round(varexp.raw[dimX],1),"%)")) + ylab(paste0("Dim",dimY," (var: ",round(varexp.raw[dimY],1),"%)")) +
  ggtitle(paste0("PCA noob")) +
  scale_color_discrete(name = "Sex")+
  theme_minimal() #ggtheme 

grid.arrange(grobs = align_patches(list(p1,p2,p3,p4)), ncol = 2)

dev.off()


###################################################################
# Write data
###################################################################

betas <- rownames_to_column(as.data.frame(betas), var = "ID")
betas.f <- rownames_to_column(as.data.frame(betas.f), var = "ID")
#betas.f2 <- rownames_to_column(as.data.frame(betas.f2), var = "ID")

fwrite(betas,file = file.path(tables_dir, paste0(script_prefix,"_","betas_noob.txt.gz")), col.names = T, row.names = F, sep = "\t", quote = F)
fwrite(betas.f,file = file.path(tables_dir, paste0(script_prefix,"_","betas_noob_filtered.txt.gz")), col.names = T, row.names = F, sep = "\t", quote = F)
#fwrite(betas.f2,file = file.path(tables_dir, paste0(script_prefix,"_","betas_noob_filtered_gapped.txt.gz")), col.names = T, row.names = F, sep = "\t", quote = F)

# Also array annotation
fwrite(annotationEPIC,file = file.path(tables_dir, paste0(script_prefix,"_","EPICv2_annotation.txt.gz")), col.names = T, row.names = F, sep = "\t", quote = F)

# Also phenodata
fwrite(phenodata,file = file.path(tables_dir, paste0("02_phenodata.txt")), col.names = T, row.names = F, sep = "\t", quote = F)

# Save the environment
save.image(paste(basedir, "/scr/02_SCRIPT_SUMMARY_FilterAndAnnot_Macroad935_Environment.RData", sep = "/"))

