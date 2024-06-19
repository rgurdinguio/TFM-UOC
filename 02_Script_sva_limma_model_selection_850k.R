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

# Download and install the required packages:
#install.packages("quadprog")
#source("http://bioconductor.org/biocLite.R")
#biocLite("limma")
#biocLite("minfi")
#biocLite("IlluminaHumanMethylationEPICmanifest")
#biocLite("IlluminaHumanMethylation450kmanifest")
#biocLite("RColorBrewer")
#biocLite("missMethyl") 
#biocLite("matrixStats")
#biocLite("Gviz") 
#biocLite("DMRcate") 
#biocLite("stringr")
#biocLite("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
#biocLite("IlluminaHumanMethylation450kanno.ilmn12.hg19")
#biocLite("GO.db")
#biocLite("FlowSorted.Blood.450k")
#biocLite("genefilter")
#biocLite("ChAMP")
#biocLite("sva")
#install.packages("FactoMineR")
#install.packages("devtools")
#install.packages("gplots")
#install.packages("cluster")
#install.packages("gridExtra")
#devtools::install_github("kassambara/factoextra")
#biocLite("ENmix")

# Load the required packages:
library("RColorBrewer")
library("pheatmap")
library("genefilter")
library(tidyr)
library(dplyr)
library(sva)
library(RColorBrewer)
library(pheatmap)
library(matrixStats)
library(minfi)
library(limma)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(DMRcate)
library(missMethyl)
library(Gviz)
library(DMRcate)

#source("~/Documents/ScriptsR/densityBeanPlot.R")

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
idatdir <- file.path(basedir, "Raw_data_macro/")

# Set up the directory for the database filtering tables (multimapping, crossreactive, SNPs...):
dataDirectory <- file.path(basedir, "Databases/Methylation_filters/")

# Set up the directory for the results:
resultsDirectory <- file.path(basedir, "PDF/")
tablesDirectory <- file.path(basedir, "Tables/")

# Load the data:
setwd(basedir)
phenoData_28_filtered <- read.delim(file = paste(tablesDirectory, "phenoData_28.txt", sep = "/"))
# phenoData_28_filtered <- read.table(file = paste(tablesDirectory, "phenoData_28.txt", sep = "/"), header = T, sep = "\t")
                                 #, stringsAsFactors = F)
# We have included a column "Group" in the database "targets" to define the groups we want to compare 
# and ordered the samples accordingly (by those groups)
# If it was not previosly done, prepare a new column "Group" and define the groups we want to compare
# phenoData_filtered$Group <- phenoData_filtered$Sample_Group 
# phenoData_filtered$Group <- gsub("Recid","R",phenoData_filtered$Group)
# phenoData_filtered$Group <- gsub("Recur","R",phenoData_filtered$Group)

# Give names to the rows of phenoData_filtered and columns of bVals and mVals tables
rownames(phenoData_28_filtered) <- phenoData_28_filtered$Condition
bVals_28 <- read.table(file = paste(tablesDirectory, "bVals_28.txt", sep = "/"), header = T, sep = "\t")
colnames(bVals_28) <- phenoData_28_filtered$Condition
mVals_28 <- read.table(file = paste(tablesDirectory, "mVals_28.txt", sep = "/"), header = T, sep = "\t")
colnames(mVals_28) <- phenoData_28_filtered$Condition


# This code will be used to select the samples that we are goint to include in the model.
# Select Progenitors and B cells:
# phenoData_filtered_2 <- phenoData_filtered[which(phenoData_filtered$Sample_Group == "B-PROGENITORS" | phenoData_filtered$Group == "MLL-AF4"),]
# bVals_filtered_2 <- bVals[,phenoData_filtered_2$Condition]
# mVals_filtered_2 <- mVals[,phenoData_filtered_2$Condition]

# Set the color palette for visual representations:
pal2 <- brewer.pal(2, "Set1")


#######################################################################################
######            INTER-SAMPLE VARIABILITY ASSESMENT AND NORMALIZATION           ######
######              ANALYSIS OF THE SOURCES OF VARIATION IN THE DATA             ######
#######################  (SCRIPT TOMADO DE github GUSTAVO)  ###########################
#######################################################################################
# Haremos gráficos de las SVAs para acotar las que se tienen en cuenta en el modelo 
# e incluir sólo las que resulten de interés (2, 3, 4...)
#
# From: https://github.com/kEYEOH/SVCONFOUND/blob/develop/vignettes/simple_use_case.Rmd
#######################################################################################

library(limma)
library(svconfound)

# A subset of our mValues dataset (5 probes x 3 samples) looks like this:
knitr::kable(mVals_28[1:5, 1:3])

# We check the phenoData_filtered characteristics:
class(phenoData_28_filtered)
dim(phenoData_28_filtered)
# A look to the first 6 rows in *phenoData*:
knitr::kable(head(phenoData_28_filtered, 6))

## SVD-based confounder analysis

# The simplest method for confounder analysis in the *svconfound* package is based
# on Singular Value Decomposition (SVD) and Principal Components Analysis (PCA), and 
# it is implemented in the *svd_analysis()* function. We can use the function
# on our data like this:

# svd_result = svd_analysis(as.data.frame(mVals), as.data.frame(phenoData_filtered))
## ERROR: Error in `contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]) : 
##        contrasts can be applied only to factors with 2 or more levels

# This error occurs when you attempt to fit a regression model using a predictor 
# variable that is either a factor or character and only has one unique value.
# We can actually use the following syntax to count the number of unique values 
# for each variable in our data frame:
# Count unique values for each variable
sapply(lapply(phenoData_28_filtered, unique), length)

# In our case, the variable "Tipo_tumoral" is considered as factor but it has only
# one value: "MAHNF". We will transform it as.character()
phenoData_28_filtered$Tipo_tumoral <- as.character(phenoData_28_filtered$Tipo_tumoral)
#phenoData_28_filtered$Comentarios <- as.character(phenoData_28_filtered$Comentarios) # NO es necesario

## If that transformation does not work, we could omit the "problematic" variables:
# phenoData_filtered2 <- phenoData_filtered[,!colnames(phenoData_filtered) == "Tipo_tumoral"]
## OR # phenoData_filtered2 <- phenoData_filtered[,!colnames(phenoData_filtered) %in% c("Sample_Plate","Tipo_tumoral")]
## Count unique values for each variable
# sapply(lapply(phenoData_filtered2, unique), length)
svd_result_28 = svd_analysis(as.data.frame(mVals_28), as.data.frame(phenoData_28_filtered))

# The *svd_analysis()* function estimates the principal components based on our 
# data, and then tests each of the components against each of the phenotype
# variables, thus producing an association score which can help the researcher to 
# determine which variables are most correlated with the greatest variations in 
# the dataset.

# We can generate a screeplot to see how much variance explains each component
# using the correspondent overloaded function:

# fig.width = 4, fig.height = 3, fig.align = 'center'
screeplot(svd_result_28)


# The darker vertical line in the plot shows the estimated number of significant
# components able to explain most of the variability present in the dataset. If 
# the researcher does not like the default plot theme, she can access the
# data used for its generation just by accessing the *variance_explained* and 
# *limit_significant_PC* fields of the resulting object.

knitr::kable(head(svd_result_28$variance_explained, 10))

svd_result_28$limit_significant_PC

pdf(file= file.path(resultsDirectory, "svd_results_28.pdf"), paper="a4r", height=21, width=28, onefile=TRUE)
par(mfrow=c(1,1))
# We can generate the main confounder plot by using the overloaded *plot()* function:
# r fig.width = 5, fig.height = 3, fig.align = 'center'}
plot(svd_result_28)

screeplot(svd_result_28)

dev.off()


# For each intersection of a principal component (PC-X) and a phenotype variable, 
# a square is shown if there is a statistically significant relationship between 
# them. The color of the square determines the significance level for their 
# relationship.

# In this case, we can see that the main source of variability (the first 
# principal component) is highly correlated with our phenotype of interest 
# (Sample_Group), which is something good and expected in this case, as this data 
# is extracted from a real dataset comparing normal and tumor tissues, where very
# high differences are expected.

# Problem is, there are other variables associated with the principal components.
# A very interesting one, from a confounder analysis point of view, is *Slide*. 
# This is a technical variable which usually affects those samples processed in 
# batches in DNA Methylation analysis pipelines. This type of confounder can 
# introduce noise in our analysis methods and should be taken into account in 
# order to obtain more accurate results.

# As in the screeplot example, the researcher can access the data used to generate
# the plot in the *significance* field.

knitr::kable(head(svd_result_28$significance, 10))


## SVA-based confounder analysis

# Another approach provided by the *svconfound* package is the one implemented in 
# the *sva_analysis()* function. It uses Surrogate Variable Analysis (SVA) in 
# order to obtain the main sources of variability in the residuals of a previously
# fitted model. Then, a similar approach to the *svd_analysis()* function is used 
# to determine the degree of association between Surrogate Variables and phenotype
# variables.

# Let's use the *sva_analysis()* function on our data:
# sva_result = sva_analysis(as.data.frame(mVals), as.data.frame(phenoData_filtered), ~ Group)

sva_result_28 = sva_analysis(as.matrix(mVals_28), as.data.frame(phenoData_28_filtered), ~ Group) 

# We have used a formula to describe our hypothesis. We are using a simple model
# containing only the *Sample_Group* variable. Our intention is to model the 
# possible sources of confounding by including the Surrogate Variables in our 
# model at a later step. We can generate a graph similar to the previous one:

pdf(file= file.path(resultsDirectory, "sva_results_28.pdf"), paper="a4r", height=21, width=28, onefile=TRUE)

# fig.width = 5, fig.height = 3, fig.align = 'center'
plot(sva_result_28)

dev.off()


# The function has estimated that the number of Surrogate Variables is 7, and the
# association graph shows that, after fitting the model, there is still a lot of 
# confounding associated to the *Slide*, *Sex* and *Age* variables. If there were
# no sources of confounding, we could go along with our simple model, but the
# presence of these confounders indicates that we could have problems if we went
# that way in this case.

# As in previous functions, the researcher can access the underlying graph
# information using the *significance* field.
knitr::kable(head(sva_result_28$significance, 10))
num_sv <- sva_result_28[["num_sv"]]

#####################################################################################
##################### SVA CORRECTION USING GUSTAVO'S SCRIPT #########################

## Fitting a model

# With all the information we have gathered, we are ready now to fit a model to 
# each of the probes, adjusting for the possible confounders we have identified
# using the previous plots. First, we should create a design matrix:
design = model.matrix(~ Group, data = phenoData_28_filtered)


# In order to correct for confounding, we are going to add the Surrogate Variables
# to our design matrix, including them as predictors in our model. 
design_sv = cbind(design, sva_result_28$surrogates)

# Now we are ready for model fitting. We are going to use *limma* for this, as it
# is a very common method in the context of microarray analyses.
fit = lmFit(mVals_28, design_sv)
efit = eBayes(fit)

# Once the model is fitted, we can obtain the p-values and effect sizes for all 
# the probes.
sig_table = topTable(
  efit, 
  coef = 'GroupR', 
  p.value = 1, 
  number = Inf
)


# With this information, we can apply whatever criteria in order to select a 
# subset of significant probes. For example, we are going to select those probes
# with an adjusted p-value below 0.38 (and an absolute effect size greater than 1.4).
sig_probes = sig_table[sig_table$adj.P.Val < 0.38 & abs(sig_table$logFC) > 0,]

##### OBTENEMOS 69 PROBES SIGNIFICATIVAS


# This gives us a total of `r nrow(sig_probes)` significant probes after adjusting
# for possible confounders. We can take a look at the 10 top rows of the list:
knitr::kable(head(sig_probes, 10))



######################   SCRIPT DE Ramón  #########################

###################################################################
#####           BATCH EFFECT CORRECTION: SVA                  #####
###################################################################

mod = model.matrix(~as.factor(Group), data=phenoData_28_filtered)
colnames(mod) <- c("Intercept", "R_vs_No")
mod0 = model.matrix(~1,data=phenoData_28_filtered)

n.sv = num.sv(as.matrix(mVals_28),mod,method="leek")
svobj = sva(as.matrix(mVals_28),mod,mod0,n.sv=n.sv)
colnames(svobj$sv) <- paste(rep("col",ncol(svobj$sv)),c(1:ncol(svobj$sv)),sep="")

modSv_allsv = cbind(mod,svobj$sv) # Aquí se puede elegir el número de svas (columnas de svobj$sv) que incluimos en el modelo
mod0Sv_allsv = cbind(mod0,svobj$sv)

modSv_7sv = cbind(mod,svobj$sv[,1:7]) # Elegimos las siete primeras SVAs (columnas de svobj$sv) en el modelo
mod0Sv_7sv = cbind(mod0,svobj$sv[,1:7])


###################################################################################
#####      PROBE-WISE DIFFERENTIAL METHYLATION ANALYSIS (SVA CORRECTION)      #####
###################################################################################

## SVA correction: considering all SVAs

# Fit the linear model with all svas:
fit_all <- lmFit(mVals_28, modSv_allsv)
R_vs_No <- factor(phenoData_28_filtered$Group)

# Create a contrast matrix for specific comparisons:
contMatrix_all <- makeContrasts(R_vs_No, levels = modSv_allsv)

# Fit the contrasts:
fit2_all <- contrasts.fit(fit_all, contMatrix_all)
fit2_all <- eBayes(fit2_all)

# Look at the numbers of Differentially Methylated CpGs at FDR < 0.05:
summary(decideTests(fit2_all))

# Get the table of results for the first contrast:
DMPs_filtered_all <- limma::topTable(fit2_all, number = Inf, coef = 1, p.value = 0.38)
DMPs_filtered_all$Methylstat <- NA
DMPs_filtered_all$Methylstat[which(DMPs_filtered_all$logFC > 0)] <- "hyper"
DMPs_filtered_all$Methylstat[which(DMPs_filtered_all$logFC < 0)] <- "hypo"


## Con este script la función num.sv determina que hay 25 SVAs. Cuando ajustamos 
## y hacemos el análisis con las 25 (_all) obtenemos p-valores (ajustados y no 
## ajustados) mayores a los del modelo de Gustavo.
## En este caso el número de DMPs filtradas con coef = 1 y p-valor ajustado de 0.38 es 0.


####################################################################################
## Other SVA model


## SVA correction: considering 7 first SVAs (based on num_sv from 
## sva_analysis with Gustavo's script)

# Fit the linear model with 7 first svas:
fit_7sva <- lmFit(mVals_28, modSv_7sv)
R_vs_No <- factor(phenoData_28_filtered$Group)

# Create a contrast matrix for specific comparisons:
contMatrix_7sva <- makeContrasts(R_vs_No, levels = modSv_7sv)

# Fit the contrasts:
fit2_7sva <- contrasts.fit(fit_7sva, contMatrix_7sva)
fit2_7sva <- eBayes(fit2_7sva)

# Look at the numbers of Differenti7svay Methylated CpGs at FDR < 0.05:
summary(decideTests(fit2_7sva))

# Get the table of results for the first contrast:
DMPs_filtered_7sva <- limma::topTable(fit2_7sva, number = Inf, coef = 1, p.value = 0.38)
DMPs_filtered_7sva$Methylstat <- NA
DMPs_filtered_7sva$Methylstat[which(DMPs_filtered_7sva$logFC > 0)] <- "hyper"
DMPs_filtered_7sva$Methylstat[which(DMPs_filtered_7sva$logFC < 0)] <- "hypo"

## Con este script la función num.sv determina que hay 24 SVAs. Cuando ajustamos 
## y hacemos el análisis con las 7 primeras SVAs (_7sva) obtenemos p-valores más 
## parecidos al modelo de Gustavo.
## En este caso el número de DMPs filtradas con coef = 1 y p-valor ajustado de 0.38 es 85.


#####################################################################################
## Other SVA model: Ramon's script forcing to 7 svas


svobj_7 = sva(as.matrix(mVals_28),mod,mod0,n.sv=7)
colnames(svobj_7$sv) <- paste(rep("col",ncol(svobj_7$sv)),c(1:ncol(svobj_7$sv)),sep="")

modSv_all7sv = cbind(mod,svobj_7$sv) # Aquí se puede elegir el número de svas (columnas de svobj$sv) que incluimos en el modelo
mod0Sv_all7sv = cbind(mod0,svobj_7$sv)

## SVA correction: considering 7 SVAs (based on sva analysis with Gustavo's script)

# Fit the linear model with 7 forced svas:
fit_all7sva <- lmFit(mVals_28, modSv_all7sv)
R_vs_No <- factor(phenoData_28_filtered$Group)

# Create a contrast matrix for specific comparisons:
contMatrix_all7sva <- makeContrasts(R_vs_No, levels = modSv_all7sv)

# Fit the contrasts:
fit2_all7sva <- contrasts.fit(fit_all7sva, contMatrix_all7sva)
fit2_all7sva <- eBayes(fit2_all7sva)

# Look at the numbers of Differenti7svay Methylated CpGs at FDR < 0.05:
summary(decideTests(fit2_all7sva))

# Get the table of results for the first contrast:
DMPs_filtered_all7sva <- limma::topTable(fit2_all7sva, number = Inf, coef = 1, p.value = 0.38)
DMPs_filtered_all7sva$Methylstat <- NA
DMPs_filtered_all7sva$Methylstat[which(DMPs_filtered_all7sva$logFC > 0)] <- "hyper"
DMPs_filtered_all7sva$Methylstat[which(DMPs_filtered_all7sva$logFC < 0)] <- "hypo"


## Con este script forzamos a que ajuste la variabilidad a 7 SVAs (n.sv=7). 
## Cuando ajustamos y hacemos el análisis con las 7 SVAs (_all7sva) obtenemos 
## exactamente los mismo p-valores que el modelo de Gustavo.
## En este caso el número de DMPs filtradas con coef = 1 y p-valor ajustado de 0.38 es 69.

### Después de comparar los distintos métodos para corregir con SVAs, dedicimos 
### utilizar el último script (que tiene exactamente los mismos resultados que
### el de Gustavo) pero, posteriormente utilizaremos el p-valor NO ajustado de 0.05.

### CONTINUAMOS CON ELLO EN EL SCRIPT 03_Script_get_DMPs.R
