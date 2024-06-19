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
library("clusterProfiler")
library("ReactomePA")
library("biomaRt")
library("enrichplot")
library("DOSE")
library("UpSetR")
library("forcats")

options("width"=110) #  using R in terminal mode on Linux, making easier the inspection of data sets

# Set up the path for the project location:
basedir <- "/media/rauldiul/Expansion/1_ROCIOyJAVI/RG_Macroadenomas"
tablesDirectory <- file.path(basedir, "Tables/")
resultsDirectory <- file.path(basedir, "PDF/")

# Set working directory:
setwd(tablesDirectory)

# Load the annotated datasets:
R_vs_No_all <- read.csv("Annotated_R_vs_No_all.txt", sep = "\t")
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

# Prepare a merged file with all the genes to use it as our "universe"
# Guardamos el conjunto de CpGs que componen nuestro "universo" (todas las CpGs de annotation_EPIC, puesto que está filtrado)
annotation_EPIC <- read.csv("annotation_EPIC.txt", sep = "\t", stringsAsFactors = F)
#rownames(annotation_EPIC) <- annotation_EPIC$Name
#annotation_EPIC$Relation_to_Island <- gsub(pattern = "N_Shelf", replacement = "Shelf", x = annotation_EPIC$Relation_to_Island)
#annotation_EPIC$Relation_to_Island <- gsub(pattern = "S_Shelf", replacement = "Shelf", x = annotation_EPIC$Relation_to_Island)
#annotation_EPIC$Relation_to_Island <- gsub(pattern = "N_Shore", replacement = "Shore", x = annotation_EPIC$Relation_to_Island)
#annotation_EPIC$Relation_to_Island <- gsub(pattern = "S_Shore", replacement = "Shore", x = annotation_EPIC$Relation_to_Island)

# Eliminimos NAs y duplicados
annotated_universe_1 <- annotation_EPIC[complete.cases(annotation_EPIC$geneId),c("Name","geneId")]
annotated_universe_2 <- annotated_universe_1[!duplicated(annotated_universe_1$geneId),] 
annotated_universe_2$geneId <- as.character(annotated_universe_2$geneId)

# Eliminamos NA y genes duplicados en el geneId de los datasets
R_vs_No_all_2 <- R_vs_No_all[complete.cases(R_vs_No_all$geneId),c("Name","geneId")] 
R_vs_No_all_3 <- R_vs_No_all_2[!duplicated(R_vs_No_all_2$geneId),] 
R_vs_No_HYPER_2 <- R_vs_No_HYPER[complete.cases(R_vs_No_HYPER$geneId),c("Name","geneId")] 
R_vs_No_HYPER_3 <- R_vs_No_HYPER_2[!duplicated(R_vs_No_HYPER_2$geneId),] 
R_vs_No_hypo_2 <- R_vs_No_hypo[complete.cases(R_vs_No_hypo$geneId),c("Name","geneId")] 
R_vs_No_hypo_3 <- R_vs_No_hypo_2[!duplicated(R_vs_No_hypo_2$geneId),] 


Recid_all_2 <- Recid_all[complete.cases(Recid_all$geneId),c("Name","geneId")] 
Recid_all_3 <- Recid_all_2[!duplicated(Recid_all_2$geneId),] 
Recid_HYPER_2 <- Recid_HYPER[complete.cases(Recid_HYPER$geneId),c("Name","geneId")] 
Recid_HYPER_3 <- Recid_HYPER_2[!duplicated(Recid_HYPER_2$geneId),] 
Recid_hypo_2 <- Recid_hypo[complete.cases(Recid_hypo$geneId),c("Name","geneId")] 
Recid_hypo_3 <- Recid_hypo_2[!duplicated(Recid_hypo_2$geneId),] 


RECUR_all_2 <- RECUR_all[complete.cases(RECUR_all$geneId),c("Name","geneId")] 
RECUR_all_3 <- RECUR_all_2[!duplicated(RECUR_all_2$geneId),] 
RECUR_HYPER_2 <- RECUR_HYPER[complete.cases(RECUR_HYPER$geneId),c("Name","geneId")] 
RECUR_HYPER_3 <- RECUR_HYPER_2[!duplicated(RECUR_HYPER_2$geneId),] 
RECUR_hypo_2 <- RECUR_hypo[complete.cases(RECUR_hypo$geneId),c("Name","geneId")] 
RECUR_hypo_3 <- RECUR_hypo_2[!duplicated(RECUR_hypo_2$geneId),] 


RECUR_parcial_all_2 <- RECUR_parcial_all[complete.cases(RECUR_parcial_all$geneId),c("Name","geneId")] 
RECUR_parcial_all_3 <- RECUR_parcial_all_2[!duplicated(RECUR_parcial_all_2$geneId),] 
RECUR_parcial_HYPER_2 <- RECUR_parcial_HYPER[complete.cases(RECUR_parcial_HYPER$geneId),c("Name","geneId")] 
RECUR_parcial_HYPER_3 <- RECUR_parcial_HYPER_2[!duplicated(RECUR_parcial_HYPER_2$geneId),] 
RECUR_parcial_hypo_2 <- RECUR_parcial_hypo[complete.cases(RECUR_parcial_hypo$geneId),c("Name","geneId")] 
RECUR_parcial_hypo_3 <- RECUR_parcial_hypo_2[!duplicated(RECUR_parcial_hypo_2$geneId),] 


# Transform the geneId in character (NOT NECESSARY)
#R_vs_No_all$geneId <- as.character(R_vs_No_all$geneId)
#R_vs_No_HYPER$geneId <- as.character(R_vs_No_HYPER$geneId)
#R_vs_No_hypo$geneId <- as.character(R_vs_No_hypo$geneId)

#Recid_all$geneId <- as.character(Recid_all$geneId)
#Recid_HYPER$geneId <- as.character(Recid_HYPER$geneId)
#Recid_hypo$geneId <- as.character(Recid_hypo$geneId)

#RECUR_all$geneId <- as.character(RECUR_all$geneId)
#RECUR_HYPER$geneId <- as.character(RECUR_HYPER$geneId)
#RECUR_hypo$geneId <- as.character(RECUR_hypo$geneId)

#RECUR_parcial_all$geneId <- as.character(RECUR_parcial_all$geneId)
#RECUR_parcial_HYPER$geneId <- as.character(RECUR_parcial_HYPER$geneId)
#RECUR_parcial_hypo$geneId <- as.character(RECUR_parcial_hypo$geneId)


# Generate lists of genes:
probesfiles_R = list(R_vs_No_all = R_vs_No_all_3$geneId,
                     R_vs_No_HYPER = R_vs_No_HYPER_3$geneId,
                     R_vs_No_hypo = R_vs_No_hypo_3$geneId)

probesfiles_Recid = list(Recid_all = Recid_all_3$geneId,
                         Recid_HYPER = Recid_HYPER_3$geneId,
                         Recid_hypo = Recid_hypo_3$geneId)

probesfiles_RECUR = list(RECUR_all = RECUR_all_3$geneId,
                         RECUR_HYPER = RECUR_HYPER_3$geneId,
                         RECUR_hypo = RECUR_hypo_3$geneId)

probesfiles_RECUR_parcial = list(RECUR_parcial_all = RECUR_parcial_all_3$geneId,
                                 RECUR_parcial_HYPER = RECUR_parcial_HYPER_3$geneId,
                                 RECUR_parcial_hypo = RECUR_parcial_hypo_3$geneId)

probesfiles_RvsNo_RECURparcial = list(R_vs_No_all = R_vs_No_all_3$geneId,
                                      R_vs_No_HYPER = R_vs_No_HYPER_3$geneId,
                                      R_vs_No_hypo = R_vs_No_hypo_3$geneId,
                                      RECUR_parcial_all = RECUR_parcial_all_3$geneId,
                                      RECUR_parcial_HYPER = RECUR_parcial_HYPER_3$geneId,
                                      RECUR_parcial_hypo = RECUR_parcial_hypo_3$geneId)

probesfiles_allComparisons = list(R_vs_No_all = R_vs_No_all_3$geneId,
                                  R_vs_No_HYPER = R_vs_No_HYPER_3$geneId,
                                  R_vs_No_hypo = R_vs_No_hypo_3$geneId,
                                  Recid_all = Recid_all_3$geneId,
                                  Recid_HYPER = Recid_HYPER_3$geneId,
                                  Recid_hypo = Recid_hypo_3$geneId,
                                  RECUR_all = RECUR_all_3$geneId,
                                  RECUR_HYPER = RECUR_HYPER_3$geneId,
                                  RECUR_hypo = RECUR_hypo_3$geneId,
                                  RECUR_parcial_all = RECUR_parcial_all_3$geneId,
                                  RECUR_parcial_HYPER = RECUR_parcial_HYPER_3$geneId,
                                  RECUR_parcial_hypo = RECUR_parcial_hypo_3$geneId)


probesfiles_allComparisons_4 = list(R_vs_No_all = R_vs_No_all_3$geneId,
                                    Recid_all = Recid_all_3$geneId,
                                    RECUR_all = RECUR_all_3$geneId,
                                    RECUR_parcial_all = RECUR_parcial_all_3$geneId)


probesfiles_allComparisons_8 = list(R_vs_No_HYPER = R_vs_No_HYPER_3$geneId,
                                    R_vs_No_hypo = R_vs_No_hypo_3$geneId,
                                    Recid_HYPER = Recid_HYPER_3$geneId,
                                    Recid_hypo = Recid_hypo_3$geneId,
                                    RECUR_HYPER = RECUR_HYPER_3$geneId,
                                    RECUR_hypo = RECUR_hypo_3$geneId,
                                    RECUR_parcial_HYPER = RECUR_parcial_HYPER_3$geneId,
                                    RECUR_parcial_hypo = RECUR_parcial_hypo_3$geneId)


probesfiles_talk = list(R_vs_No_HYPER = R_vs_No_HYPER_3$geneId,
                             R_vs_No_hypo = R_vs_No_hypo_3$geneId,
                             RECUR_parcial_HYPER = RECUR_parcial_HYPER_3$geneId,
                             RECUR_parcial_hypo = RECUR_parcial_hypo_3$geneId)


#Plot the GO analysis results

setwd(resultsDirectory)

pdf(file="GO_v5.pdf", paper = "a4r", height = 40, width = 40, onefile = TRUE)

compGO_R <- compareCluster(geneCluster   = probesfiles_R,
                           fun           = "enrichGO",
                           OrgDb         = "org.Hs.eg.db",
                           ont           = "BP",
                           qvalueCutoff  = 1,
                           pvalueCutoff  = 1,
                           pAdjustMethod = "BH",
                           universe = annotated_universe_2$geneId)
dotplot(compGO_R, showCategory = 10, title = "R_vs_No\nGO Enrichment Analysis BP clusters")

compGO_Recid <- compareCluster(geneCluster   = probesfiles_Recid,
                               fun           = "enrichGO",
                               OrgDb         = "org.Hs.eg.db",
                               ont           = "BP",
                               qvalueCutoff  = 1,
                               pvalueCutoff  = 1,
                               pAdjustMethod = "BH", 
                               universe = annotated_universe_2$geneId)
dotplot(compGO_Recid, showCategory = 10, title = "Recidiva_vs_No\nGO Enrichment Analysis BP clusters")


compGO_RECUR <- compareCluster(geneCluster = probesfiles_RECUR,
                               fun           = "enrichGO",
                               OrgDb         = "org.Hs.eg.db",
                               ont           = "BP",
                               qvalueCutoff  = 1,
                               pvalueCutoff  = 1,
                               pAdjustMethod = "BH", 
                               universe = annotated_universe_2$geneId)
dotplot(compGO_RECUR, showCategory = 10, title = "Recurrencia_vs_No\nGO Enrichment Analysis BP clusters")


compGO_RECUR_parcial <- compareCluster(geneCluster = probesfiles_RECUR_parcial,
                           fun           = "enrichGO",
                           OrgDb         = "org.Hs.eg.db",
                           ont           = "BP",
                           qvalueCutoff  = 1,
                           pvalueCutoff  = 1,
                           pAdjustMethod = "BH", 
                           universe = annotated_universe_2$geneId)
dotplot(compGO_RECUR_parcial, showCategory = 10, title = "Recurrencia_vs_No-Parcial\nGO Enrichment Analysis BP clusters")


compGO_RvsNo_RECURparcial <- compareCluster(geneCluster = probesfiles_RvsNo_RECURparcial,
                                       fun           = "enrichGO",
                                       OrgDb         = "org.Hs.eg.db",
                                       ont           = "BP",
                                       qvalueCutoff  = 1,
                                       pvalueCutoff  = 1,
                                       pAdjustMethod = "BH", 
                                       universe = annotated_universe_2$geneId)
dotplot(compGO_RvsNo_RECURparcial, showCategory = 10, title = "R_vs_No y Recurrencia_vs_No-Parcial\nGO Enrichment Analysis BP clusters")


compGO_allComparisons <- compareCluster(geneCluster = probesfiles_allComparisons,
                                            fun           = "enrichGO",
                                            OrgDb         = "org.Hs.eg.db",
                                            ont           = "BP",
                                            qvalueCutoff  = 1,
                                            pvalueCutoff  = 1,
                                            pAdjustMethod = "BH", 
                                            universe = annotated_universe_2$geneId)
dotplot(compGO_allComparisons, showCategory = 10, title = "All Comparisons\nGO Enrichment Analysis BP clusters")


compGO_allComparisons_4 <- compareCluster(geneCluster = probesfiles_allComparisons_4,
                                        fun           = "enrichGO",
                                        OrgDb         = "org.Hs.eg.db",
                                        ont           = "BP",
                                        qvalueCutoff  = 1,
                                        pvalueCutoff  = 1,
                                        pAdjustMethod = "BH", 
                                        universe = annotated_universe_2$geneId)
dotplot(compGO_allComparisons_4, showCategory = 10, title = "All Comparisons_4\nGO Enrichment Analysis BP clusters")


compGO_allComparisons_8 <- compareCluster(geneCluster = probesfiles_allComparisons_8,
                                        fun           = "enrichGO",
                                        OrgDb         = "org.Hs.eg.db",
                                        ont           = "BP",
                                        qvalueCutoff  = 1,
                                        pvalueCutoff  = 1,
                                        pAdjustMethod = "BH", 
                                        universe = annotated_universe_2$geneId)
dotplot(compGO_allComparisons_8, showCategory = 10, title = "All Comparisons_HYPER and hypo\nGO Enrichment Analysis BP clusters")

dev.off()

#######################  REPRESENTACIÓN  DE UNA  SELECCIÓN DE GO  ###############################

# Si interesa hacer un análisis de algún grupo de ontología concreto se haría lo siguiente:
# Recogemos las 30 GO más signicativas de cada lista para representarlas juntas
GO_list_R <- compGO_R@compareClusterResult[1:30,]
GO_list_Recid <- compGO_Recid@compareClusterResult[1:30,]
GO_list_RECUR <- compGO_RECUR@compareClusterResult[1:30,]
GO_list_RECUR_parcial <- compGO_RECUR_parcial@compareClusterResult[1:30,]
GO_list_compGO_RvsNo_RECURparcial <- compGO_RvsNo_RECURparcial@compareClusterResult[1:30,]
GO_list_allComparisons <- compGO_allComparisons@compareClusterResult[1:30,]
GO_list_allComparisons_4 <- compGO_allComparisons_4@compareClusterResult[1:30,]
GO_list_allComparisons_8 <- compGO_allComparisons_8@compareClusterResult[1:30,]

# Combinamos las GO de todas las comparaciones y eliminamos duplicados
GO_list <- rbind(GO_list_R, GO_list_Recid, GO_list_RECUR, GO_list_RECUR_parcial, GO_list_compGO_RvsNo_RECURparcial, GO_list_allComparisons, 
                 GO_list_allComparisons_4, GO_list_allComparisons_8) 
GO_list <- GO_list[!duplicated(GO_list$ID),c("ID","Description")]

# Guardamos esta lista en un archivo
write.table(GO_list, file = paste(tablesDirectory, "GO_list.txt", sep = "/"), sep = "\t", row.names = TRUE)
#GO_list <- read.csv("GO_list.txt", sep = "\t", header = T)



####### Si queremos ver TODAS LAS COMPARACIONES JUNTAS  ######## 

# GENERAMOS UN ARCHIVO PDF NUEVO PARA ESTOS GRÁFICOS
pdf(file="GO_list_graphs_v3.pdf", paper = "a4r", height = 40, width = 40, onefile = TRUE)

# Si queremos ver todas las comparaciones juntas: 
# Hay que tener en cuenta que la información está duplicada puesto que "_all" incluye "_HYPER" + "_hypo"
compGO_filtered <- compGO_allComparisons[which(compGO_allComparisons@compareClusterResult$ID %in% GO_list$ID)]
compGO_filtered <- compGO_filtered[which(compGO_filtered$Cluster == "R_vs_No_all" 
                                         |compGO_filtered$Cluster == "R_vs_No_HYPER"
                                         |compGO_filtered$Cluster == "R_vs_No_hypo"
                                         |compGO_filtered$Cluster == "Recid_all"
                                         |compGO_filtered$Cluster == "Recid_HYPER" 
                                         |compGO_filtered$Cluster == "Recid_hypo"
                                         |compGO_filtered$Cluster == "RECUR_all" 
                                         |compGO_filtered$Cluster == "RECUR_HYPER"
                                         |compGO_filtered$Cluster == "RECUR_hypo" 
                                         |compGO_filtered$Cluster == "RECUR_parcial_all"
                                         |compGO_filtered$Cluster == "RECUR_parcial_HYPER" 
                                         |compGO_filtered$Cluster == "RECUR_parcial_hypo"),]

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

# (OPCIONAL) Unimos el GO ID con la descripción en los comp_GO y en GO_list
# compGO_filtered$Name <- paste(compGO_filtered$ID, compGO_filtered$Description, sep = ":")
# GO_list$IDDes <- paste(GO_list$ID, GO_list$Description, sep = ":") 

# Eliminamos los qvalores no significativos, establecemos un Oddsratio máximo de "6" (y opcionalmente, ordenamos por qvalor)
compGO_filtered$Oddsratio[which(compGO_filtered$qvalue >= 0.05)] <- 0
compGO_filtered$Oddsratio[which(compGO_filtered$Oddsratio >= 6)] <- 6
#compGO_filtered <- compGO_filtered[order(compGO_filtered$log10qvalue, decreasing = T),]

# Hacemos el gráfico

ggplot(compGO_filtered, aes(Cluster, y =fct_relevel(Description, GO_list$Description))) + # fct_relevel ordena los Name en base al orden de GO_list$IDDes
  geom_point(aes(size=Oddsratio, colour = log10qvalue, fill = log10qvalue), shape = 21) + 
  scale_colour_gradient2(low = "#ffffff", high = "#39568C", midpoint = 0) + 
  scale_fill_gradient2(low = "#ffffff", high = "#39568C", midpoint = 0) + 
  scale_size_identity() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()



########     Si queremos ver OTRAS COMPARACIONES    ########## 

## GENERAMOS UN ARCHIVO PDF NUEVO PARA ESTOS GRÁFICOS
pdf(file="GO_list_graphs_v4.pdf", paper = "a4r", height = 40, width = 40, onefile = TRUE)

# Si interesa hacer un análisis de algún grupo de ontología concreto se haría lo siguiente:
#GO_list <- read.csv("GO_list.txt", sep = "\t", header = T)
compGO_filtered_1 <- compGO_allComparisons_4[which(compGO_allComparisons_4@compareClusterResult$ID %in% GO_list$ID)]
compGO_filtered_1 <- compGO_filtered_1[which(compGO_filtered_1$Cluster == "R_vs_No_all" 
                                             |compGO_filtered_1$Cluster == "Recid_all"
                                             |compGO_filtered_1$Cluster == "RECUR_all" 
                                             |compGO_filtered_1$Cluster == "RECUR_parcial_all"),]

# Calculamos el log del qvalue para mejorar su comprensión
compGO_filtered_1$log10qvalue <- -log(compGO_filtered_1$qvalue, 10)

# Calculamos el Oddsratio porque tenemos fracciones en modo texto
compGO_filtered_1 <- separate(data = compGO_filtered_1, col = "GeneRatio", into = c("a", "b"), sep = "/", remove = T)
compGO_filtered_1 <- separate(data = compGO_filtered_1, col = "BgRatio", into = c("c", "d"), sep = "/", remove = T)
compGO_filtered_1$a <- as.numeric(compGO_filtered_1$a)
compGO_filtered_1$b <- as.numeric(compGO_filtered_1$b)
compGO_filtered_1$c <- as.numeric(compGO_filtered_1$c)
compGO_filtered_1$d <- as.numeric(compGO_filtered_1$d)
compGO_filtered_1$Oddsratio <- (compGO_filtered_1$a/compGO_filtered_1$b)/(compGO_filtered_1$c/compGO_filtered_1$d)

# (OPCIONAL) Unimos el GO ID con la descripción en los comp_GO y en GO_list
# compGO_filtered_1$Name <- paste(compGO_filtered_1$ID, compGO_filtered_1$Description, sep = ":")
# GO_list$IDDes <- paste(GO_list$ID, GO_list$Description, sep = ":") 

# Eliminamos los qvalores no significativos, establecemos un Oddsratio máximo de "6" (y opcionalmente, ordenamos por qvalor)
compGO_filtered_1$Oddsratio[which(compGO_filtered_1$qvalue >= 0.05)] <- 0
compGO_filtered_1$Oddsratio[which(compGO_filtered_1$Oddsratio >= 6)] <- 6
compGO_filtered_1 <- compGO_filtered_1[order(compGO_filtered_1$log10qvalue, decreasing = T),]

# Hacemos el gráfico

ggplot(compGO_filtered_1, aes(Cluster, y =fct_relevel(Description, GO_list$Description))) + # fct_relevel ordena los Name en base al orden de GO_list$IDDes
  geom_point(aes(size=Oddsratio, colour = log10qvalue, fill = log10qvalue), shape = 21) + 
  scale_colour_gradient2(low = "#ffffff", high = "#39568C", midpoint = 0) + 
  scale_fill_gradient2(low = "#ffffff", high = "#39568C", midpoint = 0) + 
  scale_size_identity() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()





########     LA COMPARACIÓN para la charla   ########## 

## GENERAMOS UN ARCHIVO PDF NUEVO PARA ESTOS GRÁFICOS
#Plot the GO analysis results for the talk (15-06-2023)

setwd(resultsDirectory)

pdf(file="GO_talk_v1.pdf", paper = "a4r", height = 40, width = 40, onefile = TRUE)

compGO_talk <- compareCluster(geneCluster = probesfiles_talk,
                                            fun           = "enrichGO",
                                            OrgDb         = "org.Hs.eg.db",
                                            ont           = "BP",
                                            qvalueCutoff  = 1,
                                            pvalueCutoff  = 1,
                                            pAdjustMethod = "BH", 
                                            universe = annotated_universe_2$geneId)
dotplot(compGO_talk, showCategory = 10, title = "R_vs_No y Recurrencia_vs_No-Parcial HYPER e hypo\nGO Enrichment Analysis BP clusters")

dev.off()

# Si interesa hacer un análisis de algún grupo de ontología concreto se haría lo siguiente:
# Recogemos las 25 GO más signicativas de la lista para representarlas juntas
GO_list_talk <- compGO_talk@compareClusterResult[1:25,]
GO_list_talk <- GO_list_talk[!duplicated(GO_list_talk$ID),c("ID","Description")]

compGO_filtered_At <- compGO_allComparisons[which(compGO_allComparisons@compareClusterResult$ID %in% GO_list_talk$ID)]
compGO_filtered_At <- compGO_filtered_At[which(compGO_filtered_At$Cluster == "R_vs_No_HYPER"
                                         |compGO_filtered_At$Cluster == "R_vs_No_hypo"
                                         |compGO_filtered_At$Cluster == "RECUR_parcial_HYPER" 
                                         |compGO_filtered_At$Cluster == "RECUR_parcial_hypo"),]

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

# (OPCIONAL) Unimos el GO ID con la descripción en los comp_GO y en GO_list_talk
# compGO_filtered_At$Name <- paste(compGO_filtered_At$ID, compGO_filtered_At$Description, sep = ":")
# GO_list_talk$IDDes <- paste(GO_list_talk$ID, GO_list_talk$Description, sep = ":") 

# Eliminamos los qvalores no significativos, establecemos un Oddsratio máximo de "6" (y opcionalmente, ordenamos por qvalor)
compGO_filtered_At$Oddsratio[which(compGO_filtered_At$qvalue >= 0.05)] <- 0
compGO_filtered_At$Oddsratio[which(compGO_filtered_At$Oddsratio >= 6)] <- 6
compGO_filtered_At <- compGO_filtered_At[order(compGO_filtered_At$log10qvalue, decreasing = T),]

# Hacemos el gráfico

pdf(file="GO_list_talk_graphs_v2.pdf", paper = "a4r", height = 40, width = 40, onefile = TRUE)

ggplot(compGO_filtered_At, aes(Cluster, y =fct_relevel(Description, GO_list_talk$Description))) + # fct_relevel ordena los Name en base al orden de GO_list_talk$IDDes
  geom_point(aes(size=Oddsratio, colour = log10qvalue, fill = log10qvalue), shape = 21) + 
  scale_colour_gradient2(low = "#ffffff", high = "#FF0000", midpoint = 0) + 
  scale_fill_gradient2(low = "#ffffff", high = "#FF0000", midpoint = 0) + #scale_size_identity() + 
  theme_bw() + 
  theme(text = element_text(size = 16), axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()



################      ANTIGUO SCRIPT DE JAVI PARA HiC      #######################

# Si interesa hacer un análisis de algún grupo de ontología concreto se haría lo siguiente:
#list_ontologies <- read.csv("~/Epigenetica/Proyecto_Ratones_HiC/Data_Analysis/Datos_HiC_EpiLab/HiC_chicdiff_contrasts_GO_Mod_Javi/Ontology_list.csv", sep = "\t", header = F)
list_ontologies <- read.csv("~/Epigenetica/Proyecto_Ratones_HiC/Data_Analysis/Datos_HiC_EpiLab/HiC_chicdiff_Javi/7_GO/Ontology_lists/Ontology_list_mod7.csv", sep = ";", header = F)

compGO_filtered <- compGO[which(compGO@compareClusterResult$ID %in% list_ontologies$V1)]
compGO_filtered <- compGO_filtered[which(compGO_filtered$Cluster == "neurons_0.1_down_sym" 
                                         |compGO_filtered$Cluster == "neurons_0.1_up_sym"
                                         |compGO_filtered$Cluster == "merged_JC_VC_VC_VE_all_genes"
                                         ),]
compGO_filtered$log10qvalue <- -log(compGO_filtered$qvalue, 10)
compGO_filtered <- separate(data = compGO_filtered, col = "GeneRatio", into = c("a", "b"), sep = "/", remove = T)
compGO_filtered <- separate(data = compGO_filtered, col = "BgRatio", into = c("c", "d"), sep = "/", remove = T)
compGO_filtered$a <- as.numeric(compGO_filtered$a)
compGO_filtered$b <- as.numeric(compGO_filtered$b)
compGO_filtered$c <- as.numeric(compGO_filtered$c)
compGO_filtered$d <- as.numeric(compGO_filtered$d)
compGO_filtered$Oddsratio <- (compGO_filtered$a/compGO_filtered$b)/(compGO_filtered$c/compGO_filtered$d)
compGO_filtered$Name <- paste(compGO_filtered$ID, compGO_filtered$Description, sep = ":")

compGO_filtered$Oddsratio[which(compGO_filtered$qvalue >= 0.05)] <- 0
compGO_filtered$Oddsratio[which(compGO_filtered$Oddsratio >= 6)] <- 6
compGO_filtered <- compGO_filtered[order(compGO_filtered$log10qvalue, decreasing = T),]

list_ontologies$IDDes <- paste(list_ontologies$V1, list_ontologies$V2, sep = ":") 

ggplot(compGO_filtered, aes(Cluster, y =fct_relevel(Name, list_ontologies$IDDes))) + 
  geom_point(aes(size=Oddsratio, colour = log10qvalue, fill = log10qvalue), shape = 21) + 
  scale_colour_gradient2(low = "#ffffff", high = "#39568C", midpoint = 0) + 
  scale_fill_gradient2(low = "#ffffff", high = "#39568C", midpoint = 0) + 
  scale_size_identity() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

