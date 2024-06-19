#######################################################################
# PROJECT: Non-functioning Pituitary Macroadenomas
#######################################################################
## Differential Methylation with Sesame
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
SS_dir <- file.path(basedir, "SampleSheet/")

#Load the required packages:
#library(stats)
# Load the required packages:
library(sesame)
library(SummarizedExperiment)
library(dplyr) # for Goodness of fit
library(tidyr) # for Goodness of fit
library(ggplot2) # for pair-wise comparison plots
#install.packages("data.table") # for fwrite() function
library(data.table) # for fwrite() function
#library(readxl)
library(factoextra) # for PCA analysis
library(minfi)  # for plotCpg() function that let us plot the heatmap
#install.packages("pheatmap")
library(pheatmap) #for pheatmap() funtion to do heatmaps
#library(wheatmap)
#library(gprofiler2)

####################################################################### 
## LOAD BETAS AND PHENODATA
#######################################################################

bvals <- fread(file.path(tables_dir, "02_betas_noob.txt.gz"))
bvals_flt <- fread(file = file.path(tables_dir, "02_betas_noob_filtered.txt.gz"))

# Also array annotation
annotationEPIC <- fread(file.path(tables_dir, "02_EPICv2_annotation.txt.gz"))

# Also phenodata
phenodata <- fread(file.path(tables_dir, "02_phenodata.txt"))

# As factor
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

bvals <- as.data.frame(bvals)
rownames(bvals) <- bvals[,1]
bvals$ID <- NULL  # DELETE ID COLUMN (it also could be done with: bvals <- bvals[,-1])
bvals <- as.matrix(bvals)

bvals_flt <- as.data.frame(bvals_flt)
rownames(bvals_flt) <- bvals_flt[,1]
bvals_flt$ID <- NULL  # DELETE ID COLUMN (it also could be done with: bvals_flt <- bvals_flt[,-1])
bvals_flt <- as.matrix(bvals_flt)

#######################################################################
### Plot PCAs
#######################################################################

# Plot PCAs  
pdf(file=paste(results_dir, "PCA_plots_flt.pdf", sep = "/"), paper="a4r", height=21, width=28, 
    onefile=TRUE)

PCA <- prcomp(t(bvals_flt), scale. = FALSE)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

fviz_pca_ind(PCA, repel = T, geom.ind = c("point", "text"), 
             fill.ind = c(as.character(phenodata$Group)), 
             title= "Individuals PCA - Group", col.ind = "black",
             pointshape = 21, pointsize = 4, addEllipses = F, palette = "jco")

fviz_pca_ind(PCA, repel = T, geom.ind = c("point", "text"), 
             fill.ind = c(as.character(phenodata$Hospital)), 
             title= "Individuals PCA - Hospital", col.ind = "black", 
             pointshape = 21, pointsize = 4, addEllipses = F, palette = "jco")

fviz_pca_ind(PCA, repel = T, geom.ind = c("point", "text"), 
             fill.ind = c(as.character(phenodata$Sex)), 
             title= "Individuals PCA - Sex", 
             col.ind = "black", pointshape = 21, pointsize = 4, addEllipses = F, 
             palette = "jco")

fviz_pca_ind(PCA, repel = T, geom.ind = c("point", "text"), 
             fill.ind = c(as.character(phenodata$Diagnosis_Age)), 
             title = "Individuals PCA - Diagnosis_Age",
             col.ind = "black", pointshape = 21, pointsize = 4, addEllipses = F, 
             palette = "jco")

fviz_pca_ind(PCA, repel = T, geom.ind = c("point", "text"), 
             fill.ind = c(as.character(phenodata$Diagnosis_Age_Interval)), 
             title = "Individuals PCA - Diagnosis_Age_Interval",
             col.ind = "black", pointshape = 21, pointsize = 4, addEllipses = F, 
             palette = "jco")

fviz_pca_ind(PCA, repel = T, geom.ind = c("point", "text"), 
             fill.ind = c(as.character(phenodata$Histological_Subclassification)), 
             title= "Individuals PCA - Histological_Subclassification", 
             col.ind = "black", pointshape = 21, pointsize = 4, addEllipses = F, 
             palette = "jco")

fviz_pca_ind(PCA, repel = T, geom.ind = c("point", "text"), 
             fill.ind = c(as.character(phenodata$Invasion_grade_Knosp_scale)), 
             title = "Individuals PCA - Invasion_grade_Knosp_scale",
             col.ind = "black", pointshape = 21, pointsize = 4, addEllipses = F, 
             palette = "jco")

fviz_pca_ind(PCA, repel = T, geom.ind = c("point", "text"), 
             fill.ind = c(as.character(phenodata$Cavernous_sinus_invasion)), 
             title = "Individuals PCA - Cavernous_sinus_invasion",
             col.ind = "black", pointshape = 21, pointsize = 4, addEllipses = F, 
             palette = "jco")

fviz_pca_ind(PCA, repel = T, geom.ind = c("point", "text"), 
             fill.ind = c(as.character(phenodata$Resection_grade)), 
             title = "Individuals PCA - Resection_grade",
             col.ind = "black", pointshape = 21, pointsize = 4, addEllipses = F, 
             palette = "jco")

fviz_pca_ind(PCA, repel = T, geom.ind = c("point", "text"), 
             fill.ind = c(as.character(phenodata$Comorb_HTA)), 
             title= "Individuals PCA - Comorb_HTA",
             col.ind = "black", pointshape = 21, pointsize = 4, addEllipses = F, 
             palette = "jco")

fviz_pca_ind(PCA, repel = T, geom.ind = c("point", "text"), 
             fill.ind = c(as.character(phenodata$Comorb_Hypopituitarism)), 
             title= "Individuals PCA - Comorb_Hypopituitarism",
             col.ind = "black", pointshape = 21, pointsize = 4, addEllipses = F, 
             palette = "jco")

fviz_pca_ind(PCA, repel = T, geom.ind = c("point", "text"), 
             fill.ind = c(as.character(phenodata$Comorb_Hypothyroidism)), 
             title= "Individuals PCA - Comorb_Hypothyroidism",
             col.ind = "black", pointshape = 21, pointsize = 4, addEllipses = F, 
             palette = "jco")

fviz_pca_ind(PCA, repel = T, geom.ind = c("point", "text"), 
             fill.ind = c(as.character(phenodata$Smoking)), 
             title= "Individuals PCA - Smoking", col.ind = "black",
             pointshape = 21, pointsize = 4, addEllipses = F, palette = "jco")


fviz_pca_ind(PCA, repel = T, geom.ind = c("point", "text"), 
             fill.ind = c(as.character(phenodata$Deceased)), 
             title= "Individuals PCA - Deceased",
             col.ind = "black", pointshape = 21, pointsize = 4, addEllipses = F, 
             palette = "jco")

fviz_pca_ind(PCA, repel = T, geom.ind = c("point", "text"), 
             fill.ind = c(as.character(phenodata$IHQ_FT)), 
             title= "Individuals PCA - IHQ_FT",
             col.ind = "black", pointshape = 21, pointsize = 4, addEllipses = F, 
             palette = "jco")

dev.off()

#######################################################################
## Differential Methylation (after Preprocessing with Sesame)
#######################################################################

#The dataset is stored in a `SummarizedExperiment` object, which contains a data matrix combined with column-wise metadata accessible with `colData`:
se <- SummarizedExperiment(assays=list(betas=bvals_flt), colData=phenodata)
cd = as.data.frame(colData(se)); rownames(cd) = NULL
head(cd)

# CRITICAL:** If your data contains `NA`, it is required that you exclude CpGs with missing levels. For example, one cannot assess sex-specific DNA methylation for a CpG that only has non-NA measurement on one sex. Exclusion of such CpGs for differential methylation modeling can be done using the `checkLevels` function. Here, we will check this for **Group**:  
se_ok = (checkLevels(assay(se), colData(se)$Group)) #&
            #checkLevels(assay(se), colData(se)$Sex))
sum(se_ok)                      # the number of CpGs that passes
se = se[se_ok,]

# If your model include discrete contrast variables like tissue and sex, you should 
# consider explicitly turning it into a factor variable with a reference level. 
# Indicate the reference level
colData(se)$Group <- relevel(factor(colData(se)$Group), "No")
colData(se)$Recur_Relap <- relevel(factor(colData(se)$Recur_Relap), "NO")
colData(se)$sex <- relevel(factor(colData(se)$Sex), "M")
colData(se)$Invasion_grade_Knosp_scale <- relevel(factor(colData(se)$Invasion_grade_Knosp_scale), "0")
colData(se)$Cavernous_sinus_invasion <- relevel(factor(colData(se)$Cavernous_sinus_invasion), "No")
colData(se)$Resection_grade <- relevel(factor(colData(se)$Resection_grade), "Complete")
colData(se)$Comorb_HTA <- relevel(factor(colData(se)$Comorb_HTA), "No")
colData(se)$Comorb_Hypopituitarism <- relevel(factor(colData(se)$Comorb_Hypopituitarism), "No")
colData(se)$Comorb_Hypothyroidism <- relevel(factor(colData(se)$Comorb_Hypothyroidism), "No")
colData(se)$Comorb_Diabetes <- relevel(factor(colData(se)$Comorb_Diabetes), "No")
colData(se)$Comorb_other <- relevel(factor(colData(se)$Comorb_other), "No")
colData(se)$Smoking <- relevel(factor(colData(se)$Smoking), "No")

colData(se)$Deceased <- relevel(factor(colData(se)$Deceased), "No") # Only one level
colData(se)$Transplant_or_Ttm <- relevel(factor(colData(se)$Transplant_or_Ttm), "No") # Only one level

colData(se)$FSH_range <- relevel(factor(colData(se)$FSH_range), "0")
colData(se)$LH_range <- relevel(factor(colData(se)$LH_range), "0")
colData(se)$TSH_range <- relevel(factor(colData(se)$TSH_range), "0")

colData(se)$IHQ_FT <- relevel(factor(colData(se)$IHQ_FT), "NULO") # Only one level
colData(se)$Array_version <- relevel(factor(colData(se)$Array_version), "935k") # Only one level

# Then we will model DNA methylation variation treating *Group* (and *sex*) as covariates. 
# To do that we will call the `DML` function and specify the R formula `~Group` *(+ sex)*. 
# This function fits DNA methylation reading to a linear model and perform the corresponding 
# slope test and goodness-of-fit test (F-test holding out each contrast variable). See also 
# [formula](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/formula.html) for 
# how to specify lm/glm-like symbolic form for regression. All these results are returned 
# in an object of class `DMLSummary`:  

smry = DML(se, ~Group, BPPARAM = BiocParallel::SnowParam(4))   
smry   #All these results are returned in an object of class DMLSummary
# BiocParallel::SnowParam for parallelize with n cores in Windows
# You can use `mc.cores` argument to parallelize the computing.

### Test Interpretation
# The `DMLSummary` object is a list of slightly expanded `summary.lm` objects as is typically returned by `summary(lm())`. The `summaryExtractTest` function extracts some key test statistics from the `DMLSummary` object and store them in a data frame. Rows of the data frame correspond to CpGs/loci and columns contain the slopes and p-values of each variable.
test_result = summaryExtractTest(smry)
colnames(test_result) # the column names, show four groups of statistics
head(test_result)

# With the exception of the `Intercept`, there are four groups of columns, each starting with "Est_", "Pval_", "FPval_", and "Eff_" as prefix. Here are what they represent:
# **Est_\*** : The slope estimate (aka the $\beta$ coefficient, not to be confused with the DNA methylation $\beta$-value though) for continuous variable. DNA methylation difference of the current level with respect to the reference level for nominal contrast variables. Each suffix is concatenated from the contrast variable name (e.g., tissue, sex) and the level name if the contrast variable is discrete (e.g, Cecum, Esophagus, Fat).  For example, `Est_tissueFat` should be interpreted as the estimated methylation level difference of Fat compared to the reference tissue (which is `Colon`, as set above). If reference is not set, the first level in the alphabetic order is used as the reference level.  There is a special column named ``Est_`(Intercept)` ``.  It corresponds to the base-level methylation of the reference (in this case a Female Colon sample).
# **Pval_\*** : The unadjusted p-values of t-testing the slope. This represents the statistical significance of the methylation difference. For example, `Pval_GroupR` tests whether `R` is significantly different from `No` (the reference level) in DNA methylation. The ``Pval_`(Intercept)` `` tests whether the reference level is significantly different from zero.
# **FPval_\*** : The unadjusted p-value of the F-test contrasting the full model against a reduced model with the labeled contrast variable held out. Note that "Pval_" and "FPval_" are equivalent when the contrast variable is a 2-level factor, i.e., in the case of a pairwise comparison.
# **Eff_\*** : The effect size of each normial contrast variable. This is equivalent to the maximum slope subtracted by the minimum level including the reference level (0).

#Multiple-testing adjustment can be done afterwards using R's `p.adjust` function. It is integrated to the `DMR` function by default (see below).
# Calculate ajusted P-value 
test_result$FPval_Adj_Group <- p.adjust(test_result$FPval_Group, method = "fdr", 
                                        n = length(test_result$FPval_Group))

# Create a variable indicating the difference in betas between groups:
Group_No <- phenodata$Customer_M_ID[which(phenodata$Group == "No")]
Group_R <- phenodata$Customer_M_ID[which(phenodata$Group == "R")]
bvals_flt_diff <- as.data.frame(bvals_flt)
bvals_flt_diff$diff <- rowMeans(bvals_flt_diff[, Group_R]) - rowMeans(bvals_flt_diff[, Group_No])  ### Check out every time
test_result$diff <- bvals_flt_diff$diff[match(test_result$Probe_ID, rownames(bvals_flt_diff))]

# Create a variable indicating hyper or hypomethylation:
test_result$Methylstat <- NA
test_result$Methylstat[which(test_result$diff > 0)] <- "hyper"
test_result$Methylstat[which(test_result$diff < 0)] <- "hypo"

# Write table
#setwd(tables_dir)
write.table(test_result, file = paste(tables_dir, "test_result.txt", sep = "/"), 
            sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)

# Save the environment
save.image(paste(basedir, "scr/03_SCRIPT_SUMMARY_Modeling_Macroad935_v2_Environment.RData"))

#######################################################################
### Goodness of Fit
#######################################################################

# One may want to ask a question like
# > Is the CpG methylation tissue-specific?
#    rather than
# > Is the CpG more methylated in fat compared to liver?
# The first question ask about the contrast variable as a whole while the second question concerns only a specific level in the contrast variable. To answer this question, we can use an F-test contasting the full model with a reduced model with the target contrast held out. By default, this statistics is already computed in the `DML` function. The test result is recorded in the **FPval_** columns. For example, to get all CpGs that are methylated specific to sex,


test_result_RvsNo_UnAdjP005_b10 <- test_result %>% dplyr::filter(FPval_Group < 0.05, abs(diff) > 0.1) %>%
  select(Probe_ID, FPval_Group, Eff_Group, diff, Methylstat)
write.table(test_result_RvsNo_UnAdjP005_b10, file = paste(tables_dir, "test_result_RvsNo_UnAdjP005_b10.txt", sep = "/"), 
            sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)
# Here we used 0.1 as the effect size threshold. This means DNA methylation difference under 0.1 (10%) is considered not biologically meaningful. This can be a valid assumption for homogenous cell composition as most cells would be either biallelically methylated, unmethylated or monoallelically methylated. But different threshold can be used in different analysis scenarios.

test_result_RvsNo_UnAdjP005_b20 <- test_result %>% dplyr::filter(FPval_Group < 0.05, abs(diff) > 0.2) %>%
  select(Probe_ID, FPval_Group, Eff_Group, diff, Methylstat)
write.table(test_result_RvsNo_UnAdjP005_b20, file = paste(tables_dir, "test_result_RvsNo_UnAdjP005_b20.txt", sep = "/"), 
            sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)

######### Pairwise Comparison (Not necessary) ############ 
  
# From the test result, we can also ask whether the DNA methylation is different between two sexes or between two specific tissues. For example, `Est_sexMale` compares male from females. The following code creates a volcano plot.
library(ggplot2)
ggplot(test_result) + geom_point(aes(Est_GroupR, -log10(Pval_GroupR)))

#### Pairwise Comparison
test_result_GroupSig <- test_result[which(test_result$Pval_GroupR < 0.05),]
ggplot(test_result_GroupSig, aes(Est_GroupR, -log10(Pval_GroupR))) + geom_point()
# Likewise, we can ask whether DNA methylation might be different between fat and colon. We can do*  
#ggplot(test_result) + geom_point(aes(Est_tissueFat, -log10(Pval_tissueFat)))

########### Continuous Predictors (NOT OUR CASE)
# The variable tested in the `DML` function can be continuous. Suppose we are interested 
#in `Diagnosis_Age` besides `Sex`. We will call the same function but with the following formula: 
smry2 = DML(se, ~ Diagnosis_Age, BPPARAM = BiocParallel::SnowParam(4))
test_result2 = summaryExtractTest(smry2) %>% arrange(Est_Diagnosis_Age)
# Let's verify the CpGs positively associated with age.
test_result2 %>% 
  dplyr::select(Probe_ID, Est_Diagnosis_Age, Pval_Diagnosis_Age) %>% 
  tail
df = data.frame(Age = colData(se)$Diagnosis_Age,
    BetaValue = assay(se)[test_result2$Probe_ID[nrow(test_result2)],])
ggplot(df, aes(Age, BetaValue)) + geom_smooth(method="lm") + geom_point()

####################################################################
###### CLUSTERING HEATMAP
####################################################################

# We are going to represent a heatmap with the probes that  have pval < 0.05 and abs(diff) > 0.1
#RvsNo_NonAdjPval <- read.csv(paste(tables_dir, "test_result_Group_PvalNonAdj.txt", sep = "/"), 
#                             sep = "\t") #, stringsAsFactors = F)
RvsNo_UnAdjP005_b10 <- test_result_RvsNo_UnAdjP005_b10 

# También filtramos por una diferencia media de betas de 0.2
bValues_R_UnAdjP005_b10 <- as.data.frame(bvals_flt[RvsNo_UnAdjP005_b10$Probe_ID,])
DMPs_R_UnAdjP005_b10 <- RvsNo_UnAdjP005_b10[(RvsNo_UnAdjP005_b10$Probe_ID) %in% rownames(bValues_R_UnAdjP005_b10),]
RvsNo_UnAdjP005_b10 == DMPs_R_UnAdjP005_b10

annotation_RvsNo <- data.frame(group = as.character(phenodata$Group), ID = as.character(phenodata$Customer_M_ID))
annotation_RvsNo_ordered <- annotation_RvsNo[order(match(phenodata$Customer_M_ID, colnames(bValues_R_UnAdjP005_b10))),]
rownames(annotation_RvsNo_ordered) <- annotation_RvsNo_ordered$ID
annotation_RvsNo_ordered$ID <- NULL
RvsNo_ordered <- as.factor(annotation_RvsNo_ordered$group)

# Write the resulting dataframe in a CSV file:
#write.table(DMPs_R_UnAdjP05_b20, file = paste(tablesDirectory, "DMPs_R_UnAdjP05_b20.txt", sep = "/"), sep = "\t", row.names = TRUE)

# Plot the top 4 most differentiall7svay methylated CpGs and a Heatmap:
pdf(file= file.path(results_dir, "topprobes_RvsNo_UnAdjP005_b10.pdf"), paper="a4r", height=21, width=28, onefile=TRUE)
par(mfrow=c(1,2))
sapply(DMPs_R_UnAdjP005_b10$Probe_ID[1:10], 
       function(cpg){plotCpg(as.matrix(bValues_R_UnAdjP005_b10[,1:ncol(bValues_R_UnAdjP005_b10)]), 
                             cpg = cpg, pheno = RvsNo_ordered)})
bValues3_R_UnAdjP005_b10 <- bValues_R_UnAdjP005_b10[DMPs_R_UnAdjP005_b10$Probe_ID,1:ncol(bValues_R_UnAdjP005_b10)-1]

pheatmap(bValues3_R_UnAdjP005_b10, cluster_rows = T, cluster_cols = T, color = colorRampPalette(c("red","black","green"))(1000), 
         annotation = annotation_RvsNo_ordered, clustering_method = "ward.D2", border_color = NA, 
         main = "Heatmap BVals top significant probes R_vs_No \n(UnAdj.P-Value = 0.05 + mean bVal diff = 0.10)")

dev.off()

# WITH A DIFFERENCE OF 0.2 IN BETA 

# We are going to represent a heatmap with the probes that  have pval < 0.05 and abs(diff) > 0.1
#RvsNo_NonAdjPval <- read.csv(paste(tables_dir, "test_result_Group_PvalNonAdj.txt", sep = "/"), 
#                             sep = "\t") #, stringsAsFactors = F)
RvsNo_UnAdjP005_b20 <- test_result_RvsNo_UnAdjP005_b20 

# También filtramos por una diferencia media de betas de 0.2
bValues_R_UnAdjP005_b20 <- as.data.frame(bvals_flt[RvsNo_UnAdjP005_b20$Probe_ID,])
DMPs_R_UnAdjP005_b20 <- RvsNo_UnAdjP005_b20[(RvsNo_UnAdjP005_b20$Probe_ID) %in% rownames(bValues_R_UnAdjP005_b20),]
RvsNo_UnAdjP005_b20 == DMPs_R_UnAdjP005_b20

annotation_RvsNo <- data.frame(group = as.character(phenodata$Group), ID = as.character(phenodata$Customer_M_ID))
annotation_RvsNo_ordered <- annotation_RvsNo[order(match(phenodata$Customer_M_ID, colnames(bValues_R_UnAdjP005_b20))),]
rownames(annotation_RvsNo_ordered) <- annotation_RvsNo_ordered$ID
annotation_RvsNo_ordered$ID <- NULL
RvsNo_ordered <- as.factor(annotation_RvsNo_ordered$group)

# Write the resulting dataframe in a CSV file:
#write.table(DMPs_R_UnAdjP05_b20, file = paste(tablesDirectory, "DMPs_R_UnAdjP05_b20.txt", sep = "/"), sep = "\t", row.names = TRUE)

# Plot the top 4 most differentiall7svay methylated CpGs and a Heatmap:
pdf(file= file.path(results_dir, "topprobes_RvsNo_UnAdjP005_b20.pdf"), paper="a4r", height=21, width=28, onefile=TRUE)
par(mfrow=c(1,2))
sapply(DMPs_R_UnAdjP005_b20$Probe_ID[1:10], 
       function(cpg){plotCpg(as.matrix(bValues_R_UnAdjP005_b20[,1:ncol(bValues_R_UnAdjP005_b20)]), 
                             cpg = cpg, pheno = RvsNo_ordered)})
bValues3_R_UnAdjP005_b20 <- bValues_R_UnAdjP005_b20[DMPs_R_UnAdjP005_b20$Probe_ID,1:ncol(bValues_R_UnAdjP005_b20)-1]

pheatmap(bValues3_R_UnAdjP005_b20, cluster_rows = T, cluster_cols = T, color = colorRampPalette(c("red","black","green"))(1000), 
         annotation = annotation_RvsNo_ordered, clustering_method = "ward.D2", border_color = NA, 
         main = "Heatmap BVals top significant probes R_vs_No \n(UnAdj.P-Value = 0.05 + mean bVal diff = 0.20)")

dev.off()

### Save Environment
save.image(paste(basedir, "scr/03_SCRIPT_SUMMARY_sesame_Macroad935_Environment.RData", sep="/"))
