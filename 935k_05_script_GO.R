###############################################################
#################         GO analysis         #################
###############################################################

# Load the required packages:
library("genefilter")
library("cluster")
library("gridExtra")
library("ggplot2")
library("RColorBrewer")
library("reshape2")
library("dplyr")
library("tidyr")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#BiocManager::install("clusterProfiler")
library("clusterProfiler") #for compareCluster() function
library("ReactomePA")
library("biomaRt")
library("enrichplot")
library("DOSE")
library("UpSetR")
#install.packages("forcats")
library("forcats") #for fct_relevel() function


options("width"=110) #  using R in terminal mode on Linux, making easier the inspection of data sets

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

# Set up the path for the project location:
tablesDirectory <- file.path(basedir, "Tables/")
resultsDirectory <- file.path(basedir, "PDF/")

# Load the annotated datasets:
R_vs_No_all <- read.csv(paste(tables_dir, "Annotated_R_vs_No_all.txt", sep = "/"), sep = "\t")
R_vs_No_HYPER <- read.csv(paste(tables_dir,"Annotated_R_vs_No_HYPER.txt", sep = "/"), sep = "\t")
R_vs_No_hypo <- read.csv(paste(tables_dir,"Annotated_R_vs_No_hypo.txt", sep = "/"), sep = "\t")

# Prepare a merged file with all the genes to use it as our "universe"
# Guardamos el conjunto de CpGs que componen nuestro "universo" (todas las CpGs de annotation_EPIC, puesto que está filtrado)
#annotation_EPIC <- read.csv("annotationEPIC.txt", sep = "\t", stringsAsFactors = F)
annotation_EPIC <- fread(file.path(tables_dir, "05_EPICv2_annotation.txt.gz"))

# Eliminimos NAs y duplicados
annotated_universe_1 <- annotation_EPIC[complete.cases(annotation_EPIC$chipseeker_geneId),c("ID2","chipseeker_geneId")]
annotated_universe_2 <- annotated_universe_1[!duplicated(annotated_universe_1$chipseeker_geneId),] 
annotated_universe_2$chipseeker_geneId <- as.character(annotated_universe_2$chipseeker_geneId)

# Eliminamos NA y genes duplicados en el chipseeker_geneId de los datasets
R_vs_No_all_2 <- R_vs_No_all[complete.cases(R_vs_No_all$chipseeker_geneId),c("ID2","chipseeker_geneId")] 
R_vs_No_all_3 <- R_vs_No_all_2[!duplicated(R_vs_No_all_2$chipseeker_geneId),] 
R_vs_No_HYPER_2 <- R_vs_No_HYPER[complete.cases(R_vs_No_HYPER$chipseeker_geneId),c("ID2","chipseeker_geneId")] 
R_vs_No_HYPER_3 <- R_vs_No_HYPER_2[!duplicated(R_vs_No_HYPER_2$chipseeker_geneId),] 
R_vs_No_hypo_2 <- R_vs_No_hypo[complete.cases(R_vs_No_hypo$chipseeker_geneId),c("ID2","chipseeker_geneId")] 
R_vs_No_hypo_3 <- R_vs_No_hypo_2[!duplicated(R_vs_No_hypo_2$chipseeker_geneId),] 

# Generate lists of genes:
probesfiles_R = list(R_vs_No_all = R_vs_No_all_3$chipseeker_geneId,
                     R_vs_No_HYPER = R_vs_No_HYPER_3$chipseeker_geneId,
                     R_vs_No_hypo = R_vs_No_hypo_3$chipseeker_geneId)

#######################  REPRESENTACIÓN  DE GO  ###############################

#Plot the GO analysis results
pdf(file= paste(resultsDirectory,"GO_v1.pdf", sep = "/"), paper = "a4r", height = 40, width = 40, onefile = TRUE)

compGO_R <- compareCluster(geneCluster   = probesfiles_R,
                           fun           = "enrichGO",
                           OrgDb         = "org.Hs.eg.db",
                           ont           = "BP",
                           qvalueCutoff  = 1,
                           pvalueCutoff  = 1,
                           pAdjustMethod = "BH",
                           universe = annotated_universe_2$chipseeker_geneId)
dotplot(compGO_R, showCategory = 10, title = "R_vs_No HYPER and hypo\nGO Enrichment Analysis BP clusters")

dev.off()

################# OTRA REPRESENTACIÓN GO con UNA  SELECCIÓN DE GO  ############################

# Si interesa hacer un análisis de algún grupo de ontología concreto se haría lo siguiente:
# Recogemos las 25 GO más signicativas de la lista para representarlas juntas
GO_list_R2 <- compGO_R@compareClusterResult[1:25,]
GO_list_R2 <- GO_list_R2[!duplicated(GO_list_R2$ID), c("ID","Description")]

compGO_filtered_2 <- compGO_R[which(compGO_R@compareClusterResult$ID %in% GO_list_R2$ID)]
compGO_filtered_2 <- compGO_filtered_2[which(compGO_filtered_2$Cluster == "R_vs_No_HYPER"
                                             |compGO_filtered_2$Cluster == "R_vs_No_hypo"),]

# Calculamos el log del qvalue para mejorar su comprensión
compGO_filtered_2$log10qvalue <- -log(compGO_filtered_2$qvalue, 10)

# Calculamos el Oddsratio porque tenemos fracciones en modo texto
compGO_filtered_2 <- separate(data = compGO_filtered_2, col = "GeneRatio", into = c("a", "b"), sep = "/", remove = T)
compGO_filtered_2 <- separate(data = compGO_filtered_2, col = "BgRatio", into = c("c", "d"), sep = "/", remove = T)
compGO_filtered_2$a <- as.numeric(compGO_filtered_2$a)
compGO_filtered_2$b <- as.numeric(compGO_filtered_2$b)
compGO_filtered_2$c <- as.numeric(compGO_filtered_2$c)
compGO_filtered_2$d <- as.numeric(compGO_filtered_2$d)
compGO_filtered_2$Oddsratio <- (compGO_filtered_2$a/compGO_filtered_2$b)/(compGO_filtered_2$c/compGO_filtered_2$d)

# Eliminamos los qvalores no significativos, establecemos un Oddsratio máximo de "6" (y opcionalmente, ordenamos por qvalor)
compGO_filtered_2$Oddsratio[which(compGO_filtered_2$qvalue >= 0.05)] <- 0
compGO_filtered_2$Oddsratio[which(compGO_filtered_2$Oddsratio >= 6)] <- 6
compGO_filtered_2 <- compGO_filtered_2[order(compGO_filtered_2$log10qvalue, decreasing = T),]

# Hacemos el gráfico

pdf(file=paste(results_dir,"GO_list_R2_graphs_v2.pdf",sep = "/"), paper = "a4r", height = 40, width = 20, onefile = TRUE)

ggplot(compGO_filtered_2, aes(Cluster, y =fct_relevel(Description, GO_list_R2$Description))) + # fct_relevel ordena los Name en base al orden de GO_list_talk$IDDes
  geom_point(aes(size=Oddsratio, colour = log10qvalue, fill = log10qvalue), shape = 21) + 
  scale_colour_gradient2(low = "#ffffff", high = "#FF0000", midpoint = 0) + 
  scale_fill_gradient2(low = "#ffffff", high = "#FF0000", midpoint = 0) + #scale_size_identity() + 
  theme_bw() + 
  theme(text = element_text(size = 16), axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#######################  REPRESENTACIÓN  DE UNA  SELECCIÓN DE GO  ###############################

# Si interesa hacer un análisis de algún grupo de ontología concreto se haría lo siguiente:
# Recogemos las 30 GO más signicativas de cada lista para representarlas juntas
GO_list_R <- compGO_R@compareClusterResult[1:30,]

# Guardamos esta lista en un archivo
write.table(GO_list_R, file = paste(tablesDirectory, "GO_list_R.txt", sep = "/"), sep = "\t", row.names = TRUE)
#GO_list <- read.csv("GO_list.txt", sep = "\t", header = T)

####################################################################
####### Si queremos hacer otro tipo de gráficos, filtrando  ######## 
############ ESTE ES EL FORMATO QUE MEJOR QUEDA ####################

# GENERAMOS UN ARCHIVO PDF NUEVO PARA ESTOS GRÁFICOS
pdf(file=paste(results_dir,"GO_list_graphs_v3.pdf", sep = "/"), paper = "a4r", height = 40, width = 40, onefile = TRUE)

# Si queremos ver todas las comparaciones juntas: 
# Hay que tener en cuenta que la información está duplicada puesto que "_all" incluye "_HYPER" + "_hypo"
compGO_filtered <- compGO_R[which(compGO_R@compareClusterResult$ID %in% GO_list_R$ID)]
compGO_filtered <- compGO_filtered[which(compGO_filtered$Cluster == "R_vs_No_HYPER"
                                         |compGO_filtered$Cluster == "R_vs_No_hypo"),]

# Calculamos el log del qvalue para mejorar su comprensión
compGO_filtered$log10qvalue <- -log(compGO_filtered$qvalue, 10)

# Calculamos el Oddsratio porque tenemos fracciones en modo texto
compGO_filtered <- separate(data = compGO_filtered, col = "GeneRatio", into = c("a", "b"), sep = "/", remove = T)
compGO_filtered <- separate(data = compGO_filtered, col = "BgRatio", into = c("c", "d"), sep = "/", remove = T)
compGO_filtered$a <- as.numeric(compGO_filtered$a)
compGO_filtered$b <- as.numeric(compGO_filtered$b)
compGO_filtered$c <- as.numeric(compGO_filtered$c)
compGO_filtered$d <- as.numeric(compGO_filtered$d)
compGO_filtered$Oddsratio <- (compGO_filtered$a/compGO_filtered$b)/(compGO_filtered$c/compGO_filtered$d)


# Eliminamos los qvalores no significativos, establecemos un Oddsratio máximo de "6" (y opcionalmente, ordenamos por qvalor)
compGO_filtered <- compGO_filtered[order(compGO_filtered$log10qvalue, decreasing = T),]

# Hacemos el gráfico

ggplot(compGO_filtered, aes(Cluster, y =fct_relevel(Description, GO_list_R$Description))) + # fct_relevel ordena los Name en base al orden de GO_list$IDDes
  geom_point(aes(size=Oddsratio, colour = log10qvalue, fill = log10qvalue), shape = 21) + 
  scale_colour_gradient2(low = "#ffffff", high = "#FF0000", midpoint = 0) + 
  scale_fill_gradient2(low = "#ffffff", high = "#FF0000", midpoint = 0) + 
  #scale_size_identity() + 
  theme_bw() + 
  theme(text = element_text(size = 15), axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

##############################################################################
########  OTRA COMPARACIÓN  

  ## GENERAMOS UN ARCHIVO PDF NUEVO PARA ESTOS GRÁFICOS
  #Plot the GO analysis results 

probesfiles_talk = list(R_vs_No_all = R_vs_No_all_3$chipseeker_geneId,
                             R_vs_No_HYPER = R_vs_No_HYPER_3$chipseeker_geneId,
                             R_vs_No_hypo = R_vs_No_hypo_3$chipseeker_geneId)
  
    
pdf(file=paste(results_dir, "GO_talk_v1.pdf", sep = "/"), paper = "a4r", height = 40, width = 40, onefile = TRUE)

compGO_talk <- compareCluster(geneCluster = probesfiles_talk,
                                   fun           = "enrichGO",
                                   OrgDb         = "org.Hs.eg.db",
                                   ont           = "BP",
                                   qvalueCutoff  = 1,
                                   pvalueCutoff  = 1,
                                   pAdjustMethod = "BH", 
                                   universe = annotated_universe_2$chipseeker_geneId)
dotplot(compGO_talk, showCategory = 10, title = "R_vs_No HYPER e hypo\nGO Enrichment Analysis BP clusters")

dev.off()

# Si interesa hacer un análisis de algún grupo de ontología concreto se haría lo siguiente:
# Recogemos las 25 GO más signicativas de la lista para representarlas juntas
GO_list_talk <- compGO_talk@compareClusterResult[1:25,]
GO_list_talk <- GO_list_talk[!duplicated(GO_list_talk$ID),c("ID","Description")]

compGO_filtered_At <- compGO_R[which(compGO_R@compareClusterResult$ID %in% GO_list_talk$ID)]
compGO_filtered_At <- compGO_filtered_At[which(compGO_filtered_At$Cluster == "R_vs_No_HYPER"
                                               |compGO_filtered_At$Cluster == "R_vs_No_hypo"),]

# Calculamos el log del qvalue para mejorar su comprensión
compGO_filtered_At$log10qvalue <- -log(compGO_filtered_At$qvalue, 10)

# Calculamos el Oddsratio porque tenemos fracciones en modo texto
compGO_filtered_At <- separate(data = compGO_filtered_At, col = "GeneRatio", into = c("a", "b"), sep = "/", remove = T)
compGO_filtered_At <- separate(data = compGO_filtered_At, col = "BgRatio", into = c("c", "d"), sep = "/", remove = T)
compGO_filtered_At$a <- as.numeric(compGO_filtered_At$a)
compGO_filtered_At$b <- as.numeric(compGO_filtered_At$b)
compGO_filtered_At$c <- as.numeric(compGO_filtered_At$c)
compGO_filtered_At$d <- as.numeric(compGO_filtered_At$d)
compGO_filtered_At$Oddsratio <- (compGO_filtered_At$a/compGO_filtered_At$b)/(compGO_filtered_At$c/compGO_filtered_At$d)

# Eliminamos los qvalores no significativos, establecemos un Oddsratio máximo de "6" (y opcionalmente, ordenamos por qvalor)
compGO_filtered_At$Oddsratio[which(compGO_filtered_At$qvalue >= 0.05)] <- 0
compGO_filtered_At$Oddsratio[which(compGO_filtered_At$Oddsratio >= 6)] <- 6
compGO_filtered_At <- compGO_filtered_At[order(compGO_filtered_At$log10qvalue, decreasing = T),]

# Hacemos el gráfico

pdf(file="GO_list_talk_graphs_v2.pdf", paper = "a4r", height = 40, width = 20, onefile = TRUE)

ggplot(compGO_filtered_At, aes(Cluster, y =fct_relevel(Description, GO_list_talk$Description))) + # fct_relevel ordena los Name en base al orden de GO_list_talk$IDDes
  geom_point(aes(size=Oddsratio, colour = log10qvalue, fill = log10qvalue), shape = 21) + 
  scale_colour_gradient2(low = "#ffffff", high = "#FF0000", midpoint = 0) + 
  scale_fill_gradient2(low = "#ffffff", high = "#FF0000", midpoint = 0) + #scale_size_identity() + 
  theme_bw() + 
  theme(text = element_text(size = 16), axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
