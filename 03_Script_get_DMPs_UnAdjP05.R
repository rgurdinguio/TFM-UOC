### PARTIMOS DEL SCRIPT 02_Script_sva_limma_model_selection.R
### Después de comparar los distintos métodos para corregir con 
### SVAs, dedicimos utilizar el último script de Ramón (que 
### tiene exactamente los mismos resultados que el de Gustavo). 

### Vamos a tomar un p-valor NO ajustado de 0.05 (a partir del 
### apartado "GET DMPs") y diff de bVals > 0.20 para mirar las 
### distintas comparaciones de grupos.

###################################################################
#####                     PREPARING DATA                      #####
###################################################################

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

options("width"=110) #  using R in terminal mode on Linux, making easier the inspection of data sets


# Set up the path for the project and the IDATs location:
basedir <- "/media/rauldiul/Expansion/1_ROCIOyJAVI/RG_Macroadenomas"
idatdir <- file.path(basedir, "Raw_data_macro/")

# Set up the directory for the database filtering tables (multimapping, crossreactive, SNPs...):
dataDirectory <- file.path(basedir, "Databases/Methylation_filters/")

# Set up the directory for the results:
resultsDirectory <- file.path(basedir, "PDF/")
tablesDirectory <- file.path(basedir, "Tables/")

# Load the data:
setwd(basedir)
phenoData_28_filtered <- read.delim("/media/rauldiul/Expansion/1_ROCIOyJAVI/RG_Macroadenomas/Tables/phenoData_28.txt")
sapply(lapply(phenoData_28_filtered, unique), length)
# Como las variables "Tipo_tumoral" y "Comentarios" son consideradas 
# factores de 1 nivel, tenemos que modificarlas:
phenoData_28_filtered$Tipo_tumoral <- as.character(phenoData_28_filtered$Tipo_tumoral)
phenoData_28_filtered$Comentarios <- as.character(phenoData_28_filtered$Comentarios)

# Give names to the rows of phenoData_filtered and columns of bVals and mVals tables
rownames(phenoData_28_filtered) <- phenoData_28_filtered$Condition
bVals_28 <- read.table(file = paste(tablesDirectory, "bVals_28.txt", sep = "/"), header = T, sep = "\t")
colnames(bVals_28) <- phenoData_28_filtered$Condition
mVals_28 <- read.table(file = paste(tablesDirectory, "mVals_28.txt", sep = "/"), header = T, sep = "\t")
colnames(mVals_28) <- phenoData_28_filtered$Condition

# Set the color palette for visual representations:
pal2 <- brewer.pal(2, "Set1")



##################################################################
###  Vamos a usar el script de Ramón para corregir por SVAs    ###
###   (forzando a num_sv = 7, de acuerdo a los resultados      ###
###        obtenidos con la función sva_analysis() )           ###
##################################################################



###########################################################################################
##################     COMPARISON-1:  No  VS  Recidiva+RECURRENCIA      ###################
###########################################################################################

mod = model.matrix(~as.factor(Group), data=phenoData_28_filtered)
colnames(mod) <- c("Intercept", "R_vs_No")
mod0 = model.matrix(~1,data=phenoData_28_filtered)

# Si quisiéramos calcular el número TOTAL de posibles SVAs, podemos hacer:
# n.sv = num.sv(as.matrix(mVals_28),mod,method="leek")

# Para analizar el mínimo número de SVAs significativas, hacemos:
library(svconfound)
sva_result_28 = sva_analysis(as.matrix(mVals_28), as.data.frame(phenoData_28_filtered), ~ Group) 

# Plot the SVAs (fig.width = 5, fig.height = 3, fig.align = 'center')
plot(sva_result_28)

# The function has estimated that the number of Surrogate Variables is 7.
num_sv_28 <- sva_result_28[["num_sv"]]

# We can access the underlying graph information using the *significance* field.
knitr::kable(head(sva_result_28$significance, 10))


#########################################################################
##########  SVA CORRECTION: considering the calculated SVAs  ############
#########################################################################

## SVA correction considering the SVAs calculated by sva_analysis() 
## (based on num_sv from sva_analysis as in Gustavo's script)
## We use the script from Ramón, forcing to the SVAs calculated in sva_analysis()
svobj_7 = sva(as.matrix(mVals_28),mod,mod0,n.sv=num_sv_28)
colnames(svobj_7$sv) <- paste(rep("col",ncol(svobj_7$sv)),c(1:ncol(svobj_7$sv)),sep="")

modSv_all7sv = cbind(mod,svobj_7$sv) # Aquí se puede elegir el número de svas (columnas de svobj$sv) que incluimos en el modelo
mod0Sv_all7sv = cbind(mod0,svobj_7$sv)

# Fit the linear model with the SVAs:
fit_all7sva <- lmFit(mVals_28, modSv_all7sv)
R_vs_No <- factor(phenoData_28_filtered$Group)

# Create a contrast matrix for specific comparisons:
contMatrix_all7sva <- makeContrasts(R_vs_No, levels = modSv_all7sv)

# Fit the contrasts:
fit2_all7sva <- contrasts.fit(fit_all7sva, contMatrix_all7sva)
fit2_all7sva <- eBayes(fit2_all7sva)

# Look at the numbers of Differentially Methylated CpGs at FDR < 0.05:
summary(decideTests(fit2_all7sva))


############################################################
############            GET DMPs               #############
############################################################

# Get the table of results for the first contrast

# Prirmero pasamos todos los valores a un objeto (data frame) 
# y filtramos por un p-valor NO ajustado de 0.05
DMPs_R <- limma::topTable(fit2_all7sva, number = Inf, coef = 1)
DMPs_R_UnAdjP05 <- DMPs_R[(DMPs_R$P.Value < 0.05),]
DMPs_R_UnAdjP05$Methylstat <- NA
DMPs_R_UnAdjP05$Methylstat[which(DMPs_R_UnAdjP05$logFC > 0)] <- "hyper"
DMPs_R_UnAdjP05$Methylstat[which(DMPs_R_UnAdjP05$logFC < 0)] <- "hypo"
# También filtramos por una diferencia media de betas de 0.2
bValues2_R_UnAdjP05 <- as.data.frame(bVals_28[rownames(DMPs_R_UnAdjP05),])
bValues2_R_UnAdjP05$diff <- rowMeans(bValues2_R_UnAdjP05[,1:22]) - rowMeans(bValues2_R_UnAdjP05[,23:28])  ### Check out every time
bValues_R_UnAdjP05_b20<- bValues2_R_UnAdjP05[which(abs(bValues2_R_UnAdjP05$diff) > 0.20),]
DMPs_R_UnAdjP05_b20 <- DMPs_R_UnAdjP05[rownames(bValues_R_UnAdjP05_b20),]

# Write the resulting dataframe in a CSV file:
write.table(DMPs_R_UnAdjP05_b20, file = paste(tablesDirectory, "DMPs_R_UnAdjP05_b20.txt", sep = "/"), sep = "\t", row.names = TRUE)

# Plot the top 4 most differentiall7svay methylated CpGs and a Heatmap:
pdf(file= file.path(resultsDirectory, "topprobes_R_vs_No_UnAdjP05_b20.pdf"), paper="a4r", height=21, width=28, onefile=TRUE)
par(mfrow=c(1,2))
sapply(rownames(DMPs_R_UnAdjP05_b20)[1:10], 
       function(cpg){plotCpg(as.matrix(bValues_R_UnAdjP05_b20[,1:ncol(bValues_R_UnAdjP05_b20)-1]), 
                             cpg = cpg, pheno = R_vs_No)})
bValues3_R_UnAdjP05_b20 <- bValues_R_UnAdjP05_b20[rownames(DMPs_R_UnAdjP05_b20),1:ncol(bValues_R_UnAdjP05_b20)-1]
annotation_R_vs_No <- data.frame(group = as.character(phenoData_28_filtered$Group))
rownames(annotation_R_vs_No) <- colnames(bValues_R_UnAdjP05_b20[,1:ncol(bValues_R_UnAdjP05_b20)-1])
pheatmap(bValues3_R_UnAdjP05_b20, cluster_rows = T, cluster_cols = T, color = colorRampPalette(c("red","black","green"))(1000), 
         annotation = annotation_R_vs_No, clustering_method = "ward.D2", border_color = NA, 
         main = "Heatmap BVals top significant probes R_vs_No \n(UnAdj.P-Value = 0.05 + mean bVal diff = 0.20)")

dev.off()

#(brewer.pal(n = 7, name = "Reds"))

###########################################################################################
####################      COMPARISON-2:   No VS Recidiva (R1)         #####################
###########################################################################################


###################################################################
#####                       READING DATA                      #####
###################################################################

# This code is used to select the samples that we are goint to include in the model.
# Select samples "No" and "Recid" in the variable "Sample_Group":
phenoData_28_filtered_R1 <- phenoData_28_filtered[which(phenoData_28_filtered$Sample_Group == "No" | 
                                                          phenoData_28_filtered$Sample_Group == "Recid"),]
bVals_28_filtered_R1 <- bVals_28[,colnames(bVals_28) %in% phenoData_28_filtered_R1$Condition]
mVals_28_filtered_R1 <- mVals_28[,colnames(mVals_28) %in% phenoData_28_filtered_R1$Condition]

## SVA model: Ramon's script forcing to the calculated SVAs
mod_R1 = model.matrix(~as.factor(Group), data=phenoData_28_filtered_R1)
colnames(mod_R1) <- c("Intercept", "R1_vs_No")
mod0_R1 = model.matrix(~1,data=phenoData_28_filtered_R1)

# Si quisiéramos calcular el número TOTAL de posibles SVAs, podemos hacer:
# n.sv = num.sv(as.matrix(mVals_28),mod,method="leek")

# Para analizar el mínimo número de SVAs significativas, hacemos:
#library(svconfound)
sva_result_R1 = sva_analysis(as.matrix(mVals_28_filtered_R1), as.data.frame(phenoData_28_filtered_R1), ~ Group) 

# Plot the SVAs (fig.width = 5, fig.height = 3, fig.align = 'center')
plot(sva_result_R1)

# The function has estimated that the number of Surrogate Variables is 6.
num_sv_R1 <- sva_result_R1[["num_sv"]]

# We can access the underlying graph information using the *significance* field.
knitr::kable(head(sva_result_R1$significance, 10))


#########################################################################
##########  SVA CORRECTION: considering the calculated SVAs  ############
#########################################################################

## SVA correction considering the SVAs calculated by sva_analysis() 
## (based on num_sv from sva_analysis as in Gustavo's script)
## We use the script from Ramón, forcing to the SVAs calculated in sva_analysis()
svobj_R1 = sva(as.matrix(mVals_28_filtered_R1),mod_R1,mod0_R1,n.sv=num_sv_R1)
colnames(svobj_R1$sv) <- paste(rep("col",ncol(svobj_R1$sv)),c(1:ncol(svobj_R1$sv)),sep="")

modSv_R1 = cbind(mod_R1,svobj_R1$sv) # Aquí se puede elegir el número de svas (columnas de svobj$sv) que incluimos en el modelo
mod0Sv_R1 = cbind(mod0_R1,svobj_R1$sv)

# Fit the linear model with 7 SVAs:
fit_R1 <- lmFit(mVals_28_filtered_R1, modSv_R1)
R1_vs_No <- factor(phenoData_28_filtered_R1$Group)

# Create a contrast matrix for specific comparisons:
contMatrix_R1 <- makeContrasts(R1_vs_No, levels = modSv_R1)

# Fit the contrasts:
fit2_R1 <- contrasts.fit(fit_R1, contMatrix_R1)
fit2_R1 <- eBayes(fit2_R1)

# Look at the numbers of Differenti7svay Methylated CpGs at FDR < 0.05:
summary(decideTests(fit2_R1))


############################################################
############            GET DMPs               #############
############################################################

# Get the table of results for the first contrast

# Prirmero pasamos todos los valores a un objeto (data frame) 
# y filtramos por un p-valor NO ajustado de 0.05
DMPs_R1 <- limma::topTable(fit2_R1, number = Inf, coef = 1)
DMPs_R1_UnAdjP05 <- DMPs_R1[(DMPs_R1$P.Value < 0.05),]
DMPs_R1_UnAdjP05$Methylstat <- NA
DMPs_R1_UnAdjP05$Methylstat[which(DMPs_R1_UnAdjP05$logFC > 0)] <- "hyper"
DMPs_R1_UnAdjP05$Methylstat[which(DMPs_R1_UnAdjP05$logFC < 0)] <- "hypo"
# También filtramos por una diferencia media de betas de 0.2
bValues2_R1_UnAdjP05 <- as.data.frame(bVals_28_filtered_R1[rownames(DMPs_R1_UnAdjP05),])
bValues2_R1_UnAdjP05$diff <- rowMeans(bValues2_R1_UnAdjP05[,1:22]) - rowMeans(bValues2_R1_UnAdjP05[,23:24])  ### Check out every time
bValues_R1_UnAdjP05_b20 <- bValues2_R1_UnAdjP05[which(abs(bValues2_R1_UnAdjP05$diff) > 0.20),]
DMPs_R1_UnAdjP05_b20 <- DMPs_R1_UnAdjP05[rownames(bValues_R1_UnAdjP05_b20),]

# Write the resulting dataframe in a CSV file:
write.table(DMPs_R1_UnAdjP05_b20, file = paste(tablesDirectory, "DMPs_R1_UnAdjP05_b20.txt", sep = "/"), sep = "\t", row.names = TRUE)

# Plot the top 4 most differentially methylated CpGs and a Heatmap:
pdf(file= file.path(resultsDirectory, "topprobes_R1_vs_No_UnAdjP05_b20.pdf"), paper="a4r", height=21, width=28, onefile=TRUE)
par(mfrow=c(1,2))
sapply(rownames(DMPs_R1_UnAdjP05_b20)[1:10], 
       function(cpg){plotCpg(as.matrix(bValues_R1_UnAdjP05_b20[,1:ncol(bValues_R1_UnAdjP05_b20)-1]), 
                             cpg = cpg, pheno = R1_vs_No)})
bValues3_R1_UnAdjP05_b20 <- bValues_R1_UnAdjP05_b20[rownames(DMPs_R1_UnAdjP05_b20),1:ncol(bValues_R1_UnAdjP05_b20)-1]
annotation_R1_vs_No <- data.frame(group = as.character(phenoData_28_filtered_R1$Group))
rownames(annotation_R1_vs_No) <- colnames(bValues_R1_UnAdjP05_b20[,1:ncol(bValues_R1_UnAdjP05_b20)-1])
pheatmap(bValues3_R1_UnAdjP05_b20, cluster_rows = T, cluster_cols = T, color = colorRampPalette(c("red","black","green"))(1000), 
         annotation = annotation_R_vs_No, clustering_method = "ward.D2", border_color = NA, 
         main = "Heatmap BVals top significant probes Recid_vs_No \n(UnAdj.P-Value = 0.05 + mean bVal diff = 0.20)")

dev.off()


###########################################################################################
####################      COMPARISON-3:   No VS RECURRENCIA (R2)      #####################
###########################################################################################


###################################################################
#####                       READING DATA                      #####
###################################################################


# This code is used to select the samples that we are goint to include in the model.
# Select samples "No" and "Recur" in the variable "Sample_Group":
phenoData_28_filtered_R2 <- phenoData_28_filtered[which(phenoData_28_filtered$Sample_Group == "No" | 
                                                          phenoData_28_filtered$Sample_Group == "Recur"),]
bVals_28_filtered_R2 <- bVals_28[,colnames(bVals_28) %in% phenoData_28_filtered_R2$Condition]
mVals_28_filtered_R2 <- mVals_28[,colnames(mVals_28) %in% phenoData_28_filtered_R2$Condition]

## SVA model: Ramon's script forcing to the calculated svas
mod_R2 = model.matrix(~as.factor(Group), data=phenoData_28_filtered_R2)
colnames(mod_R2) <- c("Intercept", "R2_vs_No")
mod0_R2 = model.matrix(~1,data=phenoData_28_filtered_R2)

# Si quisiéramos ca lcular el número TOTAL de posibles SVAs, podemos hacer:
# n.sv = num.sv(as.matrix(mVals_28),mod,method="leek")

# Para analizar el mínimo número de SVAs significativas, hacemos:
#library(svconfound)
sva_result_R2 = sva_analysis(as.matrix(mVals_28_filtered_R2), as.data.frame(phenoData_28_filtered_R2), ~ Group) 
# Da error pero indica que el Número de variables surrogadas significativas es 7.

# Plot the SVAs (fig.width = 5, fig.height = 3, fig.align = 'center')
# plot(sva_result_R2)

# The function has estimated that the number of Surrogate Variables is 7.
num_sv_R2 <- sva_result_28[["num_sv"]]

# We can access the underlying graph information using the *significance* field.
#knitr::kable(head(sva_result_R2$significance, 10))


#########################################################################
##########  SVA CORRECTION: considering the calculated SVAs  ############
#########################################################################

## SVA correction considering the SVAs calculated by sva_analysis() 
## (based on num_sv from sva_analysis as in Gustavo's script)
## We use the script from Ramón, forcing to the SVAs calculated in sva_analysis()
svobj_R2 = sva(as.matrix(mVals_28_filtered_R2),mod_R2,mod0_R2,n.sv=num_sv_R2)
colnames(svobj_R2$sv) <- paste(rep("col",ncol(svobj_R2$sv)),c(1:ncol(svobj_R2$sv)),sep="")

modSv_R2 = cbind(mod_R2,svobj_R2$sv) # Aquí se puede elegir el número de svas (columnas de svobj$sv) que incluimos en el modelo
mod0Sv_R2 = cbind(mod0_R2,svobj_R2$sv)

# Fit the linear model with 7 SVAs:
fit_R2 <- lmFit(mVals_28_filtered_R2, modSv_R2)
R2_vs_No <- factor(phenoData_28_filtered_R2$Group)

# Create a contrast matrix for specific comparisons:
contMatrix_R2 <- makeContrasts(R2_vs_No, levels = modSv_R2)

# Fit the contrasts:
fit2_R2 <- contrasts.fit(fit_R2, contMatrix_R2)
fit2_R2 <- eBayes(fit2_R2)

# Look at the numbers of Differenti7svay Methylated CpGs at FDR < 0.05:
summary(decideTests(fit2_R2))


############################################################
############            GET DMPs               #############
############################################################

# Get the table of results for the first contrast
# Prirmero pasamos todos los valores a un objeto (data frame) 
# y filtramos por un p-valor NO ajustado de 0.05
DMPs_R2 <- limma::topTable(fit2_R2, number = Inf, coef = 1)
DMPs_R2_UnAdjP05 <- DMPs_R2[(DMPs_R2$P.Value < 0.05),]
DMPs_R2_UnAdjP05$Methylstat <- NA
DMPs_R2_UnAdjP05$Methylstat[which(DMPs_R2_UnAdjP05$logFC > 0)] <- "hyper"
DMPs_R2_UnAdjP05$Methylstat[which(DMPs_R2_UnAdjP05$logFC < 0)] <- "hypo"
# También filtramos por una diferencia media de betas de 0.2
bValues2_R2_UnAdjP05 <- as.data.frame(bVals_28_filtered_R2[rownames(DMPs_R2_UnAdjP05),])
bValues2_R2_UnAdjP05$diff <- rowMeans(bValues2_R2_UnAdjP05[,1:22]) - rowMeans(bValues2_R2_UnAdjP05[,23:26])  ### Check out every time
bValues_R2_UnAdjP05_b20 <- bValues2_R2_UnAdjP05[which(abs(bValues2_R2_UnAdjP05$diff) > 0.20),]
DMPs_R2_UnAdjP05_b20 <- DMPs_R2_UnAdjP05[rownames(bValues_R2_UnAdjP05_b20),]

# Write the resulting dataframe in a CSV file:
write.table(DMPs_R2_UnAdjP05_b20, file = paste(tablesDirectory, "DMPs_R2_UnAdjP05_b20.txt", sep = "/"), sep = "\t", row.names = TRUE)

# Plot the top 4 most differentially methylated CpGs:
pdf(file= file.path(resultsDirectory, "topprobes_R2_vs_No_UnAdjP05_b20.pdf"), paper="a4r", height=21, width=28, onefile=TRUE)
par(mfrow=c(1,2))
sapply(rownames(DMPs_R2_UnAdjP05_b20)[1:10], 
       function(cpg){plotCpg(as.matrix(bValues_R2_UnAdjP05_b20[,1:ncol(bValues_R2_UnAdjP05_b20)-1]), 
                             cpg = cpg, pheno = R2_vs_No)})
bValues3_R2_UnAdjP05_b20 <- bValues_R2_UnAdjP05_b20[rownames(DMPs_R2_UnAdjP05_b20),1:ncol(bValues_R2_UnAdjP05_b20)-1]
annotation_R2_vs_No <- data.frame(group = as.character(phenoData_28_filtered_R2$Group))
rownames(annotation_R2_vs_No) <- colnames(bValues_R2_UnAdjP05_b20[,1:ncol(bValues_R2_UnAdjP05_b20)-1])
pheatmap(bValues3_R2_UnAdjP05_b20, cluster_rows = T, cluster_cols = T, color = colorRampPalette(c("red","black","green"))(1000), 
         annotation = annotation_R_vs_No, clustering_method = "ward.D2", border_color = NA, 
         main = "Heatmap BVals top significant probes Recur_vs_No \n(UnAdj.P-Value = 0.05 + mean bVal diff = 0.20)")

dev.off()


###########################################################################################
###############      COMPARISON-4:   No-Parcial VS RECURRENCIA (R2)      ##################
###########################################################################################


###################################################################
#####                       READING DATA                      #####
###################################################################


# This code is used to select the samples that we are goint to include in the model.
# Select samples "No - Parcial" and "Recur" in the variable "Sample_Group":
phenoData_28_filtered_R2_parcial <- phenoData_28_filtered[which((phenoData_28_filtered$Sample_Group == "No" & 
                                                           phenoData_28_filtered$Grado_de_reseccion == "Parcial") | 
                                                          (phenoData_28_filtered$Sample_Group == "Recur" & 
                                                             phenoData_28_filtered$Grado_de_reseccion == "Parcial")),]
bVals_28_filtered_R2_parcial <- bVals_28[,colnames(bVals_28) %in% phenoData_28_filtered_R2_parcial$Condition]
mVals_28_filtered_R2_parcial <- mVals_28[,colnames(mVals_28) %in% phenoData_28_filtered_R2_parcial$Condition]

## SVA model: Ramon's script forcing to 7 svas
mod_R2_parcial = model.matrix(~as.factor(Group), data=phenoData_28_filtered_R2_parcial)
colnames(mod_R2_parcial) <- c("Intercept", "R2p_vs_No")
mod0_R2_parcial = model.matrix(~1,data=phenoData_28_filtered_R2_parcial)

# Si quisiéramos calcular el número TOTAL de posibles SVAs, podemos hacer:
# n.sv = num.sv(as.matrix(mVals_28),mod,method="leek")

# Para analizar el mínimo número de SVAs significativas, hacemos:
#library(svconfound)
sva_result_R2p = sva_analysis(as.matrix(mVals_28_filtered_R2_parcial), as.data.frame(phenoData_28_filtered_R2_parcial), ~ Group) 
# Da error pero indica que el Número de variables surrogadas significativas es 5.

# Plot the SVAs (fig.width = 5, fig.height = 3, fig.align = 'center')
#plot(sva_result_R2p)

# The function has estimated that the number of Surrogate Variables is 5.
num_sv_R2p <- sva_result_28[["num_sv"]]-2

# We can access the underlying graph information using the *significance* field.
#knitr::kable(head(sva_result_R2p$significance, 10))


#########################################################################
##########  SVA CORRECTION: considering the calculated SVAs  ############
#########################################################################

## SVA correction considering the SVAs calculated by sva_analysis() 
## (based on num_sv from sva_analysis as in Gustavo's script)
## We use the script from Ramón, forcing to the SVAs calculated in sva_analysis()
svobj_R2_parcial = sva(as.matrix(mVals_28_filtered_R2_parcial),mod_R2_parcial,mod0_R2_parcial,n.sv=num_sv_R2p)
colnames(svobj_R2_parcial$sv) <- paste(rep("col",ncol(svobj_R2_parcial$sv)),c(1:ncol(svobj_R2_parcial$sv)),sep="")

modSv_R2_parcial = cbind(mod_R2_parcial,svobj_R2_parcial$sv) # Aquí se puede elegir el número de svas (columnas de svobj$sv) que incluimos en el modelo
mod0Sv_R2_parcial = cbind(mod0_R2_parcial,svobj_R2_parcial$sv)

# Fit the linear model with 7 SVAs:
fit_R2_parcial <- lmFit(mVals_28_filtered_R2_parcial, modSv_R2_parcial)
R2p_vs_No <- factor(phenoData_28_filtered_R2_parcial$Group)

# Create a contrast matrix for specific comparisons:
contMatrix_R2_parcial <- makeContrasts(R2p_vs_No, levels = modSv_R2_parcial)

# Fit the contrasts:
fit2_R2_parcial <- contrasts.fit(fit_R2_parcial, contMatrix_R2_parcial)
fit2_R2_parcial <- eBayes(fit2_R2_parcial)

# Look at the numbers of Differenti7svay Methylated CpGs at FDR < 0.05:
summary(decideTests(fit2_R2_parcial))


############################################################
############            GET DMPs               #############
############################################################

# Get the table of results for the first contrast
# Prirmero pasamos todos los valores a un objeto (data frame) 
# y filtramos por un p-valor NO ajustado de 0.05
DMPs_R2_parcial <- limma::topTable(fit2_R2_parcial, number = Inf, coef = 1)
DMPs_R2_parcial_UnAdjP05 <- DMPs_R2_parcial[(DMPs_R2_parcial$P.Value < 0.05),]
DMPs_R2_parcial_UnAdjP05$Methylstat <- NA
DMPs_R2_parcial_UnAdjP05$Methylstat[which(DMPs_R2_parcial_UnAdjP05$logFC > 0)] <- "hyper"
DMPs_R2_parcial_UnAdjP05$Methylstat[which(DMPs_R2_parcial_UnAdjP05$logFC < 0)] <- "hypo"
# También filtramos por una diferencia media de betas de 0.2
bValues2_R2_parcial_UnAdjP05 <- as.data.frame(bVals_28_filtered_R2_parcial[rownames(DMPs_R2_parcial_UnAdjP05),])
bValues2_R2_parcial_UnAdjP05$diff <- rowMeans(bValues2_R2_parcial_UnAdjP05[,1:12]) - rowMeans(bValues2_R2_parcial_UnAdjP05[,13:16])  ### Check out every time
bValues_R2_parcial_UnAdjP05_b20 <- bValues2_R2_parcial_UnAdjP05[which(abs(bValues2_R2_parcial_UnAdjP05$diff) > 0.20),]
DMPs_R2_parcial_UnAdjP05_b20 <- DMPs_R2_parcial_UnAdjP05[rownames(bValues_R2_parcial_UnAdjP05_b20),]

# Write the resulting dataframe in a CSV file:
write.table(DMPs_R2_parcial_UnAdjP05_b20, file = paste(tablesDirectory, "DMPs_R2_parcial_UnAdjP05_b20.txt", sep = "/"), sep = "\t", row.names = TRUE)

# Plot the top 4 most differentially methylated CpGs:
pdf(file= file.path(resultsDirectory, "topprobes_R2_parcial_vs_No_UnAdjP05_b20.pdf"), paper="a4r", height=21, width=28, onefile=TRUE)
par(mfrow=c(1,2))
sapply(rownames(DMPs_R2_parcial_UnAdjP05_b20)[1:10], 
       function(cpg){plotCpg(as.matrix(bValues_R2_parcial_UnAdjP05_b20[,1:ncol(bValues_R2_parcial_UnAdjP05_b20)-1]), 
                             cpg = cpg, pheno = R2p_vs_No)})
bValues3_R2_parcial_UnAdjP05_b20 <- bValues_R2_parcial_UnAdjP05_b20[rownames(DMPs_R2_parcial_UnAdjP05_b20),1:ncol(bValues_R2_parcial_UnAdjP05_b20)-1]
annotation_R2_parcial_vs_No <- data.frame(group = as.character(phenoData_28_filtered_R2_parcial$Group))
rownames(annotation_R2_parcial_vs_No) <- colnames(bValues_R2_parcial_UnAdjP05_b20[,1:ncol(bValues_R2_parcial_UnAdjP05_b20)-1])
pheatmap(bValues3_R2_parcial_UnAdjP05_b20, cluster_rows = T, cluster_cols = T, color = colorRampPalette(c("red","black","green"))(1000), 
         annotation = annotation_R_vs_No, clustering_method = "ward.D2", border_color = NA, 
         main = "Heatmap BVals top significant probes Recur_vs_No_parcial\n(UnAdj.P-Value = 0.05 + mean bVal diff = 0.20)")

dev.off()



###########################################################################################
###############    COMPARISON-5:   Recidiva (R1) VS RECURRENCIA (R2)     ##################
###########################################################################################


###################################################################
#####                       READING DATA                      #####
###################################################################

# This code is used to select the samples that we are goint to include in the model.
# Select samples "Recid" and "Recur" in the variable "Sample_Group":
phenoData_28_filtered_R1_R2 <- phenoData_28_filtered[which(phenoData_28_filtered$Sample_Group == "Recid" |
                                                             phenoData_28_filtered$Sample_Group == "Recur"),]
bVals_28_filtered_R1_R2 <- bVals_28[,colnames(bVals_28) %in% phenoData_28_filtered_R1_R2$Condition]
mVals_28_filtered_R1_R2 <- mVals_28[,colnames(mVals_28) %in% phenoData_28_filtered_R1_R2$Condition]

## SVA model: Ramon's script forcing to the svas
mod_R1_R2 = model.matrix(~as.factor(Recurrencia), data=phenoData_28_filtered_R1_R2)
colnames(mod_R1_R2) <- c("Intercept", "R1_vs_R2")
mod0_R1_R2 = model.matrix(~1,data=phenoData_28_filtered_R1_R2)

# Si quisiéramos calcular el número TOTAL de posibles SVAs, podemos hacer:
# n.sv = num.sv(as.matrix(mVals_28),mod,method="leek")

# Para analizar el mínimo número de SVAs significativas, hacemos:
#library(svconfound)
sva_result_R1_R2 = sva_analysis(as.matrix(mVals_28_filtered_R1_R2), as.data.frame(phenoData_28_filtered_R1_R2), ~ Recurrencia) 
# Da error pero indica que el Número de variables surrogadas significativas es 2.

# Plot the SVAs (fig.width = 5, fig.height = 3, fig.align = 'center')
#plot(sva_result_R1_R2)

# The function has estimated that the number of Surrogate Variables is 2.
num_sv_R1_R2 <- sva_result_28[["num_sv"]]-5

# We can access the underlying graph information using the *significance* field.
#knitr::kable(head(sva_result_R1_R2$significance, 10))


#########################################################################
##########  SVA CORRECTION: considering the calculated SVAs  ############
#########################################################################

## SVA correction considering the SVAs calculated by sva_analysis() 
## (based on num_sv from sva_analysis as in Gustavo's script)
## We use the script from Ramón, forcing to the SVAs calculated in sva_analysis()
svobj_R1_R2 = sva(as.matrix(mVals_28_filtered_R1_R2),mod_R1_R2,mod0_R1_R2,n.sv=num_sv_R1_R2)
colnames(svobj_R1_R2$sv) <- paste(rep("col",ncol(svobj_R1_R2$sv)),c(1:ncol(svobj_R1_R2$sv)),sep="")

modSv_R1_R2 = cbind(mod_R1_R2,svobj_R1_R2$sv) # Aquí se puede elegir el número de svas (columnas de svobj$sv) que incluimos en el modelo
mod0Sv_R1_R2 = cbind(mod0_R1_R2,svobj_R1_R2$sv)

# Fit the linear model with 7 SVAs:
fit_R1_R2 <- lmFit(mVals_28_filtered_R1_R2, modSv_R1_R2)
R2p_vs_No <- factor(phenoData_28_filtered_R1_R2$Group)

# Create a contrast matrix for specific comparisons:
contMatrix_R1_R2 <- makeContrasts(R2p_vs_No, levels = modSv_R1_R2)

# Fit the contrasts:
fit2_R1_R2 <- contrasts.fit(fit_R1_R2, contMatrix_R1_R2)
fit2_R1_R2 <- eBayes(fit2_R1_R2)

# Look at the numbers of Differenti7svay Methylated CpGs at FDR < 0.05:
summary(decideTests(fit2_R1_R2))


############################################################
############            GET DMPs  (NO)         #############
############################################################

# Get the table of results for the first contrast
# Prirmero pasamos todos los valores a un objeto (data frame) 
# y filtramos por un p-valor NO ajustado de 0.05
DMPs_R1_R2 <- limma::topTable(fit2_R1_R2, number = Inf, coef = 1)
DMPs_R1_R2_UnAdjP05 <- DMPs_R1_R2[(DMPs_R1_R2$P.Value < 0.05),]
DMPs_R1_R2_UnAdjP05$Methylstat <- NA
DMPs_R1_R2_UnAdjP05$Methylstat[which(DMPs_R1_R2_UnAdjP05$logFC > 0)] <- "hyper"
DMPs_R1_R2_UnAdjP05$Methylstat[which(DMPs_R1_R2_UnAdjP05$logFC < 0)] <- "hypo"
# También filtramos por una diferencia media de betas de 0.2
bValues2_R1_R2_UnAdjP05 <- as.data.frame(bVals_28_filtered_R1_R2[rownames(DMPs_R1_R2_UnAdjP05),])
bValues2_R1_R2_UnAdjP05$diff <- rowMeans(bValues2_R1_R2_UnAdjP05[,1:12]) - rowMeans(bValues2_R1_R2_UnAdjP05[,13:16])  ### Check out every time
bValues_R1_R2_UnAdjP05_b20 <- bValues2_R1_R2_UnAdjP05[which(abs(bValues2_R1_R2_UnAdjP05$diff) > 0.20),]
DMPs_R1_R2_UnAdjP05_b20 <- DMPs_R1_R2_UnAdjP05[rownames(bValues_R1_R2_UnAdjP05_b20),]

# Write the resulting dataframe in a CSV file:
write.table(DMPs_R1_R2_UnAdjP05_b20, file = paste(tablesDirectory, "DMPs_R1_R2_UnAdjP05_b20.txt", sep = "/"), sep = "\t", row.names = TRUE)

# Plot the top 4 most differentially methylated CpGs:
pdf(file= file.path(resultsDirectory, "topprobes_R1_R2_vs_No_UnAdjP05_b20.pdf"), paper="a4r", height=21, width=28, onefile=TRUE)
par(mfrow=c(1,2))
sapply(rownames(DMPs_R1_R2_UnAdjP05_b20)[1:10], 
       function(cpg){plotCpg(as.matrix(bValues_R1_R2_UnAdjP05_b20[,1:ncol(bValues_R1_R2_UnAdjP05_b20)-1]), 
                             cpg = cpg, pheno = R2p_vs_No)})
bValues3_R1_R2_UnAdjP05_b20 <- bValues_R1_R2_UnAdjP05_b20[rownames(DMPs_R1_R2_UnAdjP05_b20),1:ncol(bValues_R1_R2_UnAdjP05_b20)-1]
annotation_R1_R2_vs_No <- data.frame(group = as.character(phenoData_28_filtered_R1_R2$Group))
rownames(annotation_R1_R2_vs_No) <- colnames(bValues_R1_R2_UnAdjP05_b20[,1:ncol(bValues_R1_R2_UnAdjP05_b20)-1])
pheatmap(bValues3_R1_R2_UnAdjP05_b20, cluster_rows = T, cluster_cols = T, color = colorRampPalette((brewer.pal(n = 7, name = "Reds")))(100), 
         annotation = annotation_R1_R2_vs_No, clustering_method = "ward.D2", 
         main = "Heatmap BVals top significant probes Recur_vs_No_parcial \n (UnAdj.P-Value = 0.05 + mean bVal diff = 0.20)")

dev.off()



############   VENN DIAGRAMS ######################
############      EulerR     ######################

#install.packages("eulerr")

library(eulerr)

#### HYPER e Hypo

pdf(file= file.path(resultsDirectory, "Venn Diagrams_UnAdjP05.pdf"), paper="a4r", height=21, width=28, onefile=TRUE)
par(mfrow=c(3,2))

listInput_UnAdjP05_hyper <- list(R_vs_No = row.names(DMPs_R_UnAdjP05_b20[which(DMPs_R_UnAdjP05_b20$Methylstat == "hyper"),]),
                            Recid = row.names(DMPs_R1_UnAdjP05_b20[which(DMPs_R1_UnAdjP05_b20$Methylstat == "hyper"),]), 
                            #Rd_Rr = row.names(DMPs_R1_UnAdjP05_b20[which(DMPs_R1_UnAdjP05_b20$Methylstat == "hyper"),]), 
                            Recur = row.names(DMPs_R2_UnAdjP05_b20[which(DMPs_R2_UnAdjP05_b20$Methylstat == "hyper"),]))
plot(euler(listInput_UnAdjP05_hyper, shape = "ellipse"), quantities = TRUE, main = "Recid, Recur vs No - HYPERmethylation (UnAdj.P-Value = 0.05 + mean bVal diff = 0.20)")

listInput_UnAdjP05_hypo <- list(R_vs_No = row.names(DMPs_R_UnAdjP05_b20[which(DMPs_R_UnAdjP05_b20$Methylstat == "hypo"),]),
                           Recid = row.names(DMPs_R1_UnAdjP05_b20[which(DMPs_R1_UnAdjP05_b20$Methylstat == "hypo"),]), 
                           #Rd_Rc = row.names(DMPs_R1_R2_UnAdjP05_b20[which(DMPs_R1_R2_UnAdjP05_b20$Methylstat == "hypo"),]), 
                           Recur = row.names(DMPs_R2_UnAdjP05_b20[which(DMPs_R2_UnAdjP05_b20$Methylstat == "hypo"),]))
plot(euler(listInput_UnAdjP05_hypo, shape = "ellipse"), quantities = TRUE, main = "Recid, Recur vs No - Hypomethylation (UnAdj.P-Value = 0.05 + mean bVal diff = 0.20)")

dev.off()

#### R_vs_No y RECUR_parcial

pdf(file= file.path(resultsDirectory, "Venn Diagrams_UnAdjP05_v2.pdf"), paper="a4r", height=21, width=28, onefile=TRUE)
par(mfrow=c(2,1))

listInput_UnAdjP05_parcial_hyper <- list(R_vs_No = row.names(DMPs_R_UnAdjP05_b20[which(DMPs_R_UnAdjP05_b20$Methylstat == "hyper"),]),
                                 Recur_parcial = row.names(DMPs_R2_parcial_UnAdjP05_b20[which(DMPs_R2_parcial_UnAdjP05_b20$Methylstat == "hyper"),]))
plot(euler(listInput_UnAdjP05_parcial_hyper, shape = "ellipse"), quantities = TRUE, main = "R vs No / RECUR vs Parcial - HYPERmethylation (UnAdj.P-Value = 0.05 + mean bVal diff = 0.20)")

listInput_UnAdjP05_parcial_hypo <- list(R_vs_No = row.names(DMPs_R_UnAdjP05_b20[which(DMPs_R_UnAdjP05_b20$Methylstat == "hypo"),]), 
                                Recur_parcial= row.names(DMPs_R2_parcial_UnAdjP05_b20[which(DMPs_R2_parcial_UnAdjP05_b20$Methylstat == "hypo"),]))
plot(euler(listInput_UnAdjP05_parcial_hypo, shape = "ellipse"), quantities = TRUE, main = "R vs No / RECUR vs Parcial - Hypomethylation (UnAdj.P-Value = 0.05 + mean bVal diff = 0.20)")

dev.off()



#######################################################################################
## Caracterizar las CpGs de R vs No y las de Recur vs No-parcial (separar hyper/hypo)
## Caracterizar los genes correspondientes a las CpGs anteriores  (separar hyper(hypo)
