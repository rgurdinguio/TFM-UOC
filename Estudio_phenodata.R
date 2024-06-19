# Set up the path for the project:
getwd()
# If we are in "scr" folder, move up a folder by running:
setwd("..")
getwd()
basedir <- getwd()


library(ggstatsplot) #for ggbarstats() function
library(gmodels) #for CrossTable() function

# Load the data
phenoData_all <- read.csv(file = "20240617_phenodata_all_meth.csv")

#Clean and prepare phenoData_all table:
phenoData_all[,"Sample_Name_MethylArray"] <- NULL
phenoData_all[,"Condition_MethylArray"] <- NULL
phenoData_all[,"Customer_ID"] <- NULL
phenoData_all[,"Customer_name.Array_ID."] <- NULL
phenoData_all[,"Codigo_muestra"] <- NULL
phenoData_all[,"Sentrix_ID"] <- NULL
phenoData_all[,"Sentrix_Position"] <- NULL
phenoData_all[,"Otros_IDs"] <- NULL
phenoData_all[,"ID_Cris_Alicante"] <- NULL
phenoData_all[,"ID_Cris_Sevilla"] <- NULL
phenoData_all[,"Tipo_de_muestra"] <- NULL
phenoData_all[,"Concentración_ng.ul"] <- NULL

phenoData_all[,"Fecha_de_nacimiento"] <- NULL
phenoData_all[,"Fecha_de_diagnostico"] <- NULL
phenoData_all[,"Tipo_tumoral"] <- NULL
phenoData_all[,"Recurrencia"] <- NULL
phenoData_all[,"Recidiva"] <- NULL
phenoData_all[,"Recidiva_Recurrencia"] <- NULL
phenoData_all[,"Sample_Group"] <- NULL
phenoData_all[,"Fechas_de_recurrencia"] <- NULL
phenoData_all[,"Fecha_recurrencia1"] <- NULL
phenoData_all[,"Fecha_recurrencia2"] <- NULL
phenoData_all[,"Fechas_de_recidiva"] <- NULL
phenoData_all[,"Fecha_recidiva1"] <- NULL
phenoData_all[,"Fecha_recidiva2"] <- NULL

#phenoData_all[,"Ki67_percentage"] <- NULL
phenoData_all[,"Comorbilidades"] <- NULL
#phenoData_all[,"FSH"] <- NULL
#phenoData_all[,"LH"] <- NULL
phenoData_all[,"GH"] <- NULL
phenoData_all[,"GH_range"] <- NULL
phenoData_all[,"ACTH"] <- NULL
phenoData_all[,"ACTH_range"] <- NULL
phenoData_all[,"PRL"] <- NULL
phenoData_all[,"PRL_range"] <- NULL
#phenoData_all[,"TSH"] <- NULL
#phenoData_all[,"PMR"] <- NULL

phenoData_all[,"RNAseq_list2"] <- NULL
phenoData_all[,"RNAseq_list1"] <- NULL
phenoData_all[,"RNAseq_order"] <- NULL
phenoData_all[,"Tracking_Num_RNAseq"] <- NULL
phenoData_all[,"InTrack_Num_RNAseq"] <- NULL
phenoData_all[,"Total_RNA_amount"] <- NULL
phenoData_all[,"Sent_RNA_Volume"] <- NULL
phenoData_all[,"Sent_RNA_conc"] <- NULL
phenoData_all[,"GW_RNA_Qubit"] <- NULL
phenoData_all[,"GW_RNA_conc"] <- NULL
phenoData_all[,"Hospital.1"] <- NULL
phenoData_all[,"Fecha_de_nacimiento.1"] <- NULL
phenoData_all[,"Sexo.1"] <- NULL
phenoData_all[,"Edad_al_diagnostico.1"] <- NULL
phenoData_all[,"Subclasificacion_histologica.1"] <- NULL
phenoData_all[,"Recurrencia.1"] <- NULL
phenoData_all[,"Recidiva.1"] <- NULL
phenoData_all[,"Recidiva_O_Recurrencia"] <- NULL
phenoData_all[,"Recidiva_Y_Recurrencia"] <- NULL
phenoData_all[,"Fecha_de_recurrencia"] <- NULL
phenoData_all[,"Fecha_de_recidiva"] <- NULL
phenoData_all[,"Tamano_del_tumor_mm.1"] <- NULL
phenoData_all[,"Grado_de_invasion_escala_Knosp.1"] <- NULL
phenoData_all[,"Invasion_Seno_Cavernoso.1"] <- NULL
phenoData_all[,"Grado_de_reseccion.1"] <- NULL
phenoData_all[,"Ki67_porcentaje.1"] <- NULL
phenoData_all[,"Fumador.1"] <- NULL
#phenoData_all[,"Sex"] <- NULL #Tenemos Sex en otra variable
phenoData_all[,"PMR"] <- NULL
phenoData_all[,"color.1"] <- NULL
phenoData_all[,"X"] <- NULL


phenoData_all <- phenoData_all[1:65,]

# Load the data
write.table(phenoData_all, file= "20240617_phenodata_all_meth_clean.txt", 
            sep = "\t", quote = F, col.names = T, row.names = F)

#As factor
#phenoData_all$Sentrix_ID <- as.factor(phenoData_all$Sentrix_ID)
phenoData_all$Hospital <- as.factor(phenoData_all$Hospital)
phenoData_all$Sexo <- as.factor(phenoData_all$Sexo)
phenoData_all$Intervalo_de_edad <- as.factor(phenoData_all$Intervalo_de_edad)
phenoData_all$Subclasificacion_histologica <- as.factor(phenoData_all$Subclasificacion_histologica)
#phenoData_all$Recur_Relap <- as.factor(phenoData_all$Recur_Relap)
phenoData_all$Group <- as.factor(phenoData_all$Group)
phenoData_all$Grado_de_invasion_escala_Knosp <- as.factor(phenoData_all$Grado_de_invasion_escala_Knosp)
phenoData_all$Invasion_Seno_Cavernoso <- as.factor(phenoData_all$Invasion_Seno_Cavernoso)
phenoData_all$Grado_de_reseccion <- as.factor(phenoData_all$Grado_de_reseccion)
phenoData_all$Comorb_HTA <- as.factor(phenoData_all$Comorb_HTA)
phenoData_all$Comorb_Hipopituitarismo <- as.factor(phenoData_all$Comorb_Hipopituitarismo)
phenoData_all$Comorb_Hipotiroidismo <- as.factor(phenoData_all$Comorb_Hipotiroidismo)
phenoData_all$Comorb_Diabetes <- as.factor(phenoData_all$Comorb_Diabetes)
phenoData_all$Comorb_other <- as.factor(phenoData_all$Comorb_other)
phenoData_all$Fumador <- as.factor(phenoData_all$Fumador)
phenoData_all$Fallecido <- as.factor(phenoData_all$Fallecido) # Only one level
phenoData_all$Transplante_o_Tto <- as.factor(phenoData_all$Transplante_o_Tto) # Only one level
phenoData_all$FSH_range <- as.factor(phenoData_all$FSH_range)
phenoData_all$LH_range <- as.factor(phenoData_all$LH_range)
phenoData_all$TSH_range <- as.factor(phenoData_all$TSH_range)
phenoData_all$IHQ_FT <- as.factor(phenoData_all$IHQ_FT) # Only one level
phenoData_all$Array_version <- as.factor(phenoData_all$Array_version) # Only one level

phenoData_850k <- phenoData_all[which(phenoData_all$Array_version == "850k"),]
phenoData_935k <- phenoData_all[which(phenoData_all$Array_version == "935k"),]


##################################################################
##### Recurrencia/Recidiva   vs   No Recurrencia/Recidiva ########
##################################################################

#######################################################
# SEXO
# Hacemos una gráfica con los datos y los valores estadisticos de Chi-cuadrado.
ggbarstats(data = phenoData_all, x = Sexo , y = Group, bf.message = FALSE)

## Prueba exacta de Fisher
# Para saber si el Recur/Recid (sin diferenciar entre ambas) está relacionado con el Grado de resección hacemos la prueba de Fisher en la que la hipótesis nula es la igualdad de proporciones entre los dos grupos comparados: Group = NO vs Group = R.
Fisher_Sexo_all <- CrossTable(x = phenoData_all$Sexo, y = phenoData_all$Group, 
                                            fisher = TRUE, prop.r = FALSE, prop.t = FALSE, 
                                            prop.chisq = FALSE, digits = 2)
Fisher_Sexo_850k <- CrossTable(x = phenoData_850k$Sexo, y = phenoData_850k$Group, 
                                             fisher = TRUE, prop.r = FALSE, prop.t = FALSE, 
                                             prop.chisq = FALSE, digits = 2)
Fisher_Sexo_935k <- CrossTable(x = phenoData_935k$Sexo, y = phenoData_935k$Group, 
                                             fisher = TRUE, prop.r = FALSE, prop.t = FALSE, 
                                             prop.chisq = FALSE, digits = 2)

fisher_Sexo_df <- data.frame("Perc_No_H"= c(Fisher_Sexo_all$prop.col[1,1],
                                                                 Fisher_Sexo_850k$prop.col[1,1],
                                                                 Fisher_Sexo_935k$prop.col[1,1]), 
                                           "Perc_No_M"= c(Fisher_Sexo_all$prop.col[2,1],
                                                                Fisher_Sexo_850k$prop.col[2,1],
                                                                Fisher_Sexo_935k$prop.col[2,1]),
                                           "Perc_R_H"= c(Fisher_Sexo_all$prop.col[1,2],
                                                                Fisher_Sexo_850k$prop.col[1,2],
                                                                Fisher_Sexo_935k$prop.col[1,2]), 
                                           "Perc_R_M"= c(Fisher_Sexo_all$prop.col[2,2],
                                                               Fisher_Sexo_850k$prop.col[2,2],
                                                               Fisher_Sexo_935k$prop.col[2,2]),
                                           "p_value_true odds ratio is not equal to 1" = c(Fisher_Sexo_all$fisher.ts$p.value,
                                                                                   Fisher_Sexo_850k$fisher.ts$p.value,
                                                                                   Fisher_Sexo_935k$fisher.ts$p.value), 
                                           "p_value_true odds ratio is less than 1" = c(Fisher_Sexo_all$fisher.tl$p.value,
                                                                                Fisher_Sexo_850k$fisher.tl$p.value,
                                                                                Fisher_Sexo_935k$fisher.tl$p.value),
                                           "p_value_true odds ratio is greater than 1" = c(Fisher_Sexo_all$fisher.gt$p.value,
                                                                                   Fisher_Sexo_850k$fisher.gt$p.value,
                                                                                   Fisher_Sexo_935k$fisher.gt$p.value))
rownames(fisher_Sexo_df) <- c("Sexo_ALL",
                                            "Sexo_850k",
                                            "Sexo_935k")

write.table(fisher_Sexo_df, file= "20240617_data_fisher_Sexo.txt", 
            sep = "\t", quote = F, col.names = T, row.names = F)

# Si el resultado es un p-valor no significativo, NO rechazamos la hipótesis nula y NO podemos decir que los grupos son diferentes.


#######################################################
# FUMADOR
# Hacemos una gráfica con los datos y los valores estadisticos de Chi-cuadrado.
ggbarstats(data = phenoData_all, x = Fumador , y = Group, bf.message = FALSE)

## Prueba exacta de Fisher
# Para saber si el Recur/Recid (sin diferenciar entre ambas) está relacionado con el Fumador hacemos la prueba de Fisher en la que la hipótesis nula es la igualdad de proporciones entre los dos grupos comparados: Group = NO vs Group = R.
Fisher_Fumador_all <- CrossTable(x = phenoData_all$Fumador, y = phenoData_all$Group, 
                              fisher = TRUE, prop.r = FALSE, prop.t = FALSE, 
                              prop.chisq = FALSE, digits = 2)
Fisher_Fumador_850k <- CrossTable(x = phenoData_850k$Fumador, y = phenoData_850k$Group, 
                               fisher = TRUE, prop.r = FALSE, prop.t = FALSE, 
                               prop.chisq = FALSE, digits = 2)
Fisher_Fumador_935k <- CrossTable(x = phenoData_935k$Fumador, y = phenoData_935k$Group, 
                               fisher = TRUE, prop.r = FALSE, prop.t = FALSE, 
                               prop.chisq = FALSE, digits = 2)

fisher_Fumador_df <- data.frame("Perc_No_EX"= c(Fisher_Fumador_all$prop.col[1,1],
                                            Fisher_Fumador_850k$prop.col[1,1],
                                            Fisher_Fumador_935k$prop.col[1,1]), 
                             "Perc_No_No"= c(Fisher_Fumador_all$prop.col[2,1],
                                            Fisher_Fumador_850k$prop.col[2,1],
                                            Fisher_Fumador_935k$prop.col[2,1]),
                             "Perc_No_Si"= c(Fisher_Fumador_all$prop.col[3,1],
                                             Fisher_Fumador_850k$prop.col[3,1],
                                             Fisher_Fumador_935k$prop.col[3,1]),
                             "Perc_R_EX"= c(Fisher_Fumador_all$prop.col[1,2],
                                           Fisher_Fumador_850k$prop.col[1,2],
                                           Fisher_Fumador_935k$prop.col[1,2]), 
                             "Perc_R_No"= c(Fisher_Fumador_all$prop.col[2,2],
                                           Fisher_Fumador_850k$prop.col[2,2],
                                           Fisher_Fumador_935k$prop.col[2,2]),
                             "Perc_R_Si"= c(Fisher_Fumador_all$prop.col[3,2],
                                            Fisher_Fumador_850k$prop.col[3,2],
                                            Fisher_Fumador_935k$prop.col[3,2]),
                             "p_value_true odds ratio is not equal to 1" = c(Fisher_Fumador_all$fisher.ts$p.value,
                                                                             Fisher_Fumador_850k$fisher.ts$p.value,
                                                                             Fisher_Fumador_935k$fisher.ts$p.value))#, 
#                             "p_value_true odds ratio is less than 1" = c(Fisher_Fumador_all$fisher.tl$p.value,
#                                                                          Fisher_Fumador_850k$fisher.tl$p.value,
#                                                                          Fisher_Fumador_935k$fisher.tl$p.value),
#                             "p_value_true odds ratio is greater than 1" = c(Fisher_Fumador_all$fisher.gt$p.value,
#                                                                             Fisher_Fumador_850k$fisher.gt$p.value,
#                                                                             Fisher_Fumador_935k$fisher.gt$p.value))
rownames(fisher_Fumador_df) <- c("Fumador_ALL",
                              "Fumador_850k",
                              "Fumador_935k")

write.table(fisher_Fumador_df, file= "20240617_data_fisher_Fumador.txt", 
            sep = "\t", quote = F, col.names = T, row.names = F)


#######################################################
# INTERVALO DE EDAD
# Hacemos una gráfica con los datos y los valores estadisticos de Chi-cuadrado.
ggbarstats(data = phenoData_all, x = Intervalo_de_edad , y = Group, bf.message = FALSE)

## Prueba exacta de Fisher
# Para saber si el Recur/Recid (sin diferenciar entre ambas) está relacionado con el Intervalo_de_edad hacemos la prueba de Fisher en la que la hipótesis nula es la igualdad de proporciones entre los dos grupos comparados: Group = NO vs Group = R.
Fisher_Intervalo_de_edad_all <- CrossTable(x = phenoData_all$Intervalo_de_edad, y = phenoData_all$Group, 
                                 fisher = TRUE, prop.r = FALSE, prop.t = FALSE, 
                                 prop.chisq = FALSE, digits = 2)
Fisher_Intervalo_de_edad_850k <- CrossTable(x = phenoData_850k$Intervalo_de_edad, y = phenoData_850k$Group, 
                                  fisher = TRUE, prop.r = FALSE, prop.t = FALSE, 
                                  prop.chisq = FALSE, digits = 2)
Fisher_Intervalo_de_edad_935k <- CrossTable(x = phenoData_935k$Intervalo_de_edad, y = phenoData_935k$Group, 
                                  fisher = TRUE, prop.r = FALSE, prop.t = FALSE, 
                                  prop.chisq = FALSE, digits = 2)

fisher_Intervalo_de_edad_df <- data.frame("Perc_No_14_48"= c(Fisher_Intervalo_de_edad_all$prop.col[1,1],
                                                Fisher_Intervalo_de_edad_850k$prop.col[1,1],
                                                Fisher_Intervalo_de_edad_935k$prop.col[1,1]), 
                                "Perc_No_49_58"= c(Fisher_Intervalo_de_edad_all$prop.col[2,1],
                                                Fisher_Intervalo_de_edad_850k$prop.col[2,1],
                                                Fisher_Intervalo_de_edad_935k$prop.col[2,1]),
                                "Perc_No_59_68"= c(Fisher_Intervalo_de_edad_all$prop.col[3,1],
                                                Fisher_Intervalo_de_edad_850k$prop.col[3,1],
                                                Fisher_Intervalo_de_edad_935k$prop.col[3,1]),
                                "Perc_No_69_94"= c(Fisher_Intervalo_de_edad_all$prop.col[4,1],
                                                   Fisher_Intervalo_de_edad_850k$prop.col[4,1],
                                                   Fisher_Intervalo_de_edad_935k$prop.col[4,1]),
                                "Perc_R_14_48"= c(Fisher_Intervalo_de_edad_all$prop.col[1,2],
                                               Fisher_Intervalo_de_edad_850k$prop.col[1,2],
                                               Fisher_Intervalo_de_edad_935k$prop.col[1,2]), 
                                "Perc_R_49_58"= c(Fisher_Intervalo_de_edad_all$prop.col[2,2],
                                               Fisher_Intervalo_de_edad_850k$prop.col[2,2],
                                               Fisher_Intervalo_de_edad_935k$prop.col[2,2]),
                                "Perc_R_59_68"= c(Fisher_Intervalo_de_edad_all$prop.col[3,2],
                                               Fisher_Intervalo_de_edad_850k$prop.col[3,2],
                                               Fisher_Intervalo_de_edad_935k$prop.col[3,2]),
                                "Perc_R_69_94"= c(Fisher_Intervalo_de_edad_all$prop.col[4,2],
                                                   Fisher_Intervalo_de_edad_850k$prop.col[4,2],
                                                   Fisher_Intervalo_de_edad_935k$prop.col[4,2]),
                                "p_value_true odds ratio is not equal to 1" = c(Fisher_Intervalo_de_edad_all$fisher.ts$p.value,
                                                                                Fisher_Intervalo_de_edad_850k$fisher.ts$p.value,
                                                                                Fisher_Intervalo_de_edad_935k$fisher.ts$p.value))#, 
#                             "p_value_true odds ratio is less than 1" = c(Fisher_Intervalo_de_edad_all$fisher.tl$p.value,
#                                                                          Fisher_Intervalo_de_edad_850k$fisher.tl$p.value,
#                                                                          Fisher_Intervalo_de_edad_935k$fisher.tl$p.value),
#                             "p_value_true odds ratio is greater than 1" = c(Fisher_Intervalo_de_edad_all$fisher.gt$p.value,
#                                                                             Fisher_Intervalo_de_edad_850k$fisher.gt$p.value,
#                                                                             Fisher_Intervalo_de_edad_935k$fisher.gt$p.value))
rownames(fisher_Intervalo_de_edad_df) <- c("Intervalo_de_edad_ALL",
                                 "Intervalo_de_edad_850k",
                                 "Intervalo_de_edad_935k")

write.table(fisher_Intervalo_de_edad_df, file= "20240617_data_fisher_Intervalo_de_edad.txt", 
            sep = "\t", quote = F, col.names = T, row.names = F)


#######################################################
# GRADO DE RESECCIÓN
# Hacemos una gráfica con los datos y los valores estadisticos de Chi-cuadrado.
ggbarstats(data = phenoData_all, x = Grado_de_reseccion , y = Group, bf.message = FALSE)

## Prueba exacta de Fisher
# Para saber si el Recur/Recid (sin diferenciar entre ambas) está relacionado con el Grado de resección hacemos la prueba de Fisher en la que la hipótesis nula es la igualdad de proporciones entre los dos grupos comparados: Group = NO vs Group = R.
Fisher_Grado_de_reseccion_all <- CrossTable(x = phenoData_all$Grado_de_reseccion, y = phenoData_all$Group, 
                                    fisher = TRUE, prop.r = FALSE, prop.t = FALSE, 
                                    prop.chisq = FALSE, digits = 2)
Fisher_Grado_de_reseccion_850k <- CrossTable(x = phenoData_850k$Grado_de_reseccion, y = phenoData_850k$Group, 
                                         fisher = TRUE, prop.r = FALSE, prop.t = FALSE, 
                                         prop.chisq = FALSE, digits = 2)
Fisher_Grado_de_reseccion_935k <- CrossTable(x = phenoData_935k$Grado_de_reseccion, y = phenoData_935k$Group, 
                                         fisher = TRUE, prop.r = FALSE, prop.t = FALSE, 
                                         prop.chisq = FALSE, digits = 2)

fisher_Grado_de_reseccion_df <- data.frame("No_Completa_perc"= c(Fisher_Grado_de_reseccion_all$prop.col[1,1],
                                              Fisher_Grado_de_reseccion_850k$prop.col[1,1],
                                              Fisher_Grado_de_reseccion_935k$prop.col[1,1]), 
                        "No_Parcial_perc"= c(Fisher_Grado_de_reseccion_all$prop.col[2,1],
                                              Fisher_Grado_de_reseccion_850k$prop.col[2,1],
                                              Fisher_Grado_de_reseccion_935k$prop.col[2,1]),
                        "R_Completa_perc"= c(Fisher_Grado_de_reseccion_all$prop.col[1,2],
                                              Fisher_Grado_de_reseccion_850k$prop.col[1,2],
                                              Fisher_Grado_de_reseccion_935k$prop.col[1,2]), 
                        "R_Parcial_perc"= c(Fisher_Grado_de_reseccion_all$prop.col[2,2],
                                             Fisher_Grado_de_reseccion_850k$prop.col[2,2],
                                             Fisher_Grado_de_reseccion_935k$prop.col[2,2]),
                        "p_value_true odds ratio is not equal to 1" = c(Fisher_Grado_de_reseccion_all$fisher.ts$p.value,
                                                                Fisher_Grado_de_reseccion_850k$fisher.ts$p.value,
                                                                Fisher_Grado_de_reseccion_935k$fisher.ts$p.value), 
                        "p_value_true odds ratio is less than 1" = c(Fisher_Grado_de_reseccion_all$fisher.tl$p.value,
                                                             Fisher_Grado_de_reseccion_850k$fisher.tl$p.value,
                                                             Fisher_Grado_de_reseccion_935k$fisher.tl$p.value),
                        "p_value_true odds ratio is greater than 1" = c(Fisher_Grado_de_reseccion_all$fisher.gt$p.value,
                                                                Fisher_Grado_de_reseccion_850k$fisher.gt$p.value,
                                                                Fisher_Grado_de_reseccion_935k$fisher.gt$p.value))
rownames(fisher_Grado_de_reseccion_df) <- c("Grado_de_reseccion_ALL",
                                         "Grado_de_reseccion_850k",
                                         "Grado_de_reseccion_935k")

write.table(fisher_Grado_de_reseccion_df, file= "20240617_data_fisher_Grado_de_reseccion.txt", 
            sep = "\t", quote = F, col.names = T, row.names = F)


#######################################################
# SUBCLASIFICACIÓN HISTOLOGICA
# Hacemos una gráfica con los datos y los valores estadisticos de Chi-cuadrado.
ggbarstats(data = phenoData_all, x = Subclasificacion_histologica, 
           y = Group, bf.message = FALSE)

## Prueba exacta de Fisher
# Para saber si el Recur/Recid (sin diferenciar entre ambas) está relacionado con la Subclasificacion_histologica hacemos la prueba de Fisher en la que la hipótesis nula es la igualdad de proporciones entre los dos grupos comparados: Group = NO vs Group = R.
Fisher_Subclasificacion_histologica_all <- CrossTable(x = phenoData_all$Subclasificacion_histologica, y = phenoData_all$Group, 
                                            fisher = TRUE, prop.r = FALSE, prop.t = FALSE, 
                                            prop.chisq = FALSE, digits = 2)
Fisher_Subclasificacion_histologica_850k <- CrossTable(x = phenoData_850k$Subclasificacion_histologica, y = phenoData_850k$Group, 
                                             fisher = TRUE, prop.r = FALSE, prop.t = FALSE, 
                                             prop.chisq = FALSE, digits = 2)
Fisher_Subclasificacion_histologica_935k <- CrossTable(x = phenoData_935k$Subclasificacion_histologica, y = phenoData_935k$Group, 
                                             fisher = TRUE, prop.r = FALSE, prop.t = FALSE, 
                                             prop.chisq = FALSE, digits = 2)

fisher_Subclasificacion_histologica_df <- data.frame("No_CTS_perc"= c(Fisher_Subclasificacion_histologica_all$prop.col[1,1],
                                                                      Fisher_Subclasificacion_histologica_850k$prop.col[1,1],
                                                                      Fisher_Subclasificacion_histologica_935k$prop.col[1,1]), 
                                           "No_GT_perc"= c(Fisher_Subclasificacion_histologica_all$prop.col[2,1],
                                                           Fisher_Subclasificacion_histologica_850k$prop.col[2,1],
                                                           Fisher_Subclasificacion_histologica_935k$prop.col[2,1]),
                                           "No_mixto_perc"= c(Fisher_Subclasificacion_histologica_all$prop.col[3,1],
                                                              0,
                                                              Fisher_Subclasificacion_histologica_935k$prop.col[3,1]),
                                           "No_NULO_perc"= c(Fisher_Subclasificacion_histologica_all$prop.col[4,1],
                                                             Fisher_Subclasificacion_histologica_850k$prop.col[3,1],
                                                             Fisher_Subclasificacion_histologica_935k$prop.col[4,1]),
                                           "No_PHSU_perc"= c(Fisher_Subclasificacion_histologica_all$prop.col[5,1],
                                                             Fisher_Subclasificacion_histologica_850k$prop.col[4,1],
                                                             0),
                                           "R_CTS_perc"= c(Fisher_Subclasificacion_histologica_all$prop.col[1,2],
                                                            Fisher_Subclasificacion_histologica_850k$prop.col[1,2],
                                                            Fisher_Subclasificacion_histologica_935k$prop.col[1,2]), 
                                           "R_GT_perc"= c(Fisher_Subclasificacion_histologica_all$prop.col[2,2],
                                                           Fisher_Subclasificacion_histologica_850k$prop.col[2,2],
                                                           Fisher_Subclasificacion_histologica_935k$prop.col[2,2]),
                                           "R_mixto_perc"= c(Fisher_Subclasificacion_histologica_all$prop.col[3,2],
                                                              0,
                                                              Fisher_Subclasificacion_histologica_935k$prop.col[3,2]),
                                           "R_NULO_perc"= c(Fisher_Subclasificacion_histologica_all$prop.col[4,2],
                                                             Fisher_Subclasificacion_histologica_850k$prop.col[3,2],
                                                             Fisher_Subclasificacion_histologica_935k$prop.col[4,2]),
                                           "R_PHSU_perc"= c(Fisher_Subclasificacion_histologica_all$prop.col[5,2],
                                                             Fisher_Subclasificacion_histologica_850k$prop.col[4,2],
                                                             0),
                                           "p_value_true odds ratio is not equal to 1" = c(Fisher_Subclasificacion_histologica_all$fisher.ts$p.value,
                                                                                   Fisher_Subclasificacion_histologica_850k$fisher.ts$p.value,
                                                                                   Fisher_Subclasificacion_histologica_935k$fisher.ts$p.value)) #, 
#                                           "p_value_true odds ratio is less than 1" = c(Fisher_Subclasificacion_histologica_all$fisher.tl$p.value,
#                                                                                Fisher_Subclasificacion_histologica_850k$fisher.tl$p.value,
#                                                                                Fisher_Subclasificacion_histologica_935k$fisher.tl$p.value),
#                                           "p_value_true odds ratio is greater than 1" = c(Fisher_Subclasificacion_histologica_all$fisher.gt$p.value,
#                                                                                   Fisher_Subclasificacion_histologica_850k$fisher.gt$p.value,
#                                                                                   Fisher_Subclasificacion_histologica_935k$fisher.gt$p.value))
rownames(fisher_Subclasificacion_histologica_df) <- c("Subclasificacion_histologica_ALL",
                                            "Subclasificacion_histologica_850k",
                                            "Subclasificacion_histologica_935k")

write.table(fisher_Subclasificacion_histologica_df, file= "20240617_data_fisher_Subclasificacion_histologica.txt", 
            sep = "\t", quote = F, col.names = T, row.names = F)


#######################################################
# GRADO DE INVASION ESCALA KNOSP
# Hacemos una gráfica con los datos y los valores estadisticos de Chi-cuadrado.
ggbarstats(data = phenoData_all, x = Grado_de_invasion_escala_Knosp, 
           y = Group, bf.message = FALSE)

## Prueba exacta de Fisher
# Para saber si el Recur/Recid (sin diferenciar entre ambas) está relacionado con la Grado_de_invasion_escala_Knosphacemos la prueba de Fisher en la que la hipótesis nula es la igualdad de proporciones entre los dos grupos comparados: Group = NO vs Group = R.
Fisher_Grado_de_invasion_escala_Knosp_all <- CrossTable(x = phenoData_all$Grado_de_invasion_escala_Knosp, y = phenoData_all$Group, 
                                                      fisher = TRUE, prop.r = FALSE, prop.t = FALSE, 
                                                      prop.chisq = FALSE, digits = 2)
Fisher_Grado_de_invasion_escala_Knosp_850k <- CrossTable(x = phenoData_850k$Grado_de_invasion_escala_Knosp, y = phenoData_850k$Group, 
                                                       fisher = TRUE, prop.r = FALSE, prop.t = FALSE, 
                                                       prop.chisq = FALSE, digits = 2)
Fisher_Grado_de_invasion_escala_Knosp_935k <- CrossTable(x = phenoData_935k$Grado_de_invasion_escala_Knosp, y = phenoData_935k$Group, 
                                                       fisher = TRUE, prop.r = FALSE, prop.t = FALSE, 
                                                       prop.chisq = FALSE, digits = 2)

fisher_Grado_de_invasion_escala_Knosp_df <- data.frame(
  "No_0_perc"= c(Fisher_Grado_de_invasion_escala_Knosp_all$prop.col[1,1],
                 Fisher_Grado_de_invasion_escala_Knosp_850k$prop.col[1,1],
                 Fisher_Grado_de_invasion_escala_Knosp_935k$prop.col[1,1]),
  "No_1_perc"= c(Fisher_Grado_de_invasion_escala_Knosp_all$prop.col[2,1],
                 Fisher_Grado_de_invasion_escala_Knosp_850k$prop.col[2,1],
                 Fisher_Grado_de_invasion_escala_Knosp_935k$prop.col[2,1]),
  "No_2_perc"= c(Fisher_Grado_de_invasion_escala_Knosp_all$prop.col[3,1],
                 Fisher_Grado_de_invasion_escala_Knosp_850k$prop.col[3,1],
                 Fisher_Grado_de_invasion_escala_Knosp_935k$prop.col[3,1]),
  "No_3_perc"= c(Fisher_Grado_de_invasion_escala_Knosp_all$prop.col[4,1],
                 Fisher_Grado_de_invasion_escala_Knosp_850k$prop.col[4,1],
                 Fisher_Grado_de_invasion_escala_Knosp_935k$prop.col[4,1]),
  "No_4_perc"= c(Fisher_Grado_de_invasion_escala_Knosp_all$prop.col[5,1],
                 Fisher_Grado_de_invasion_escala_Knosp_850k$prop.col[5,1],
                 Fisher_Grado_de_invasion_escala_Knosp_935k$prop.col[5,1]),
  "R_0_perc"= c(Fisher_Grado_de_invasion_escala_Knosp_all$prop.col[1,2],
                Fisher_Grado_de_invasion_escala_Knosp_850k$prop.col[1,2],
                Fisher_Grado_de_invasion_escala_Knosp_935k$prop.col[1,2]), 
  "R_1_perc"= c(Fisher_Grado_de_invasion_escala_Knosp_all$prop.col[2,2],
                Fisher_Grado_de_invasion_escala_Knosp_850k$prop.col[2,2],
                Fisher_Grado_de_invasion_escala_Knosp_935k$prop.col[2,2]),
  "R_2_perc"= c(Fisher_Grado_de_invasion_escala_Knosp_all$prop.col[3,2],
                Fisher_Grado_de_invasion_escala_Knosp_850k$prop.col[3,2],
                Fisher_Grado_de_invasion_escala_Knosp_935k$prop.col[3,2]),
  "R_3_perc"= c(Fisher_Grado_de_invasion_escala_Knosp_all$prop.col[4,2],
                Fisher_Grado_de_invasion_escala_Knosp_850k$prop.col[4,2],
                Fisher_Grado_de_invasion_escala_Knosp_935k$prop.col[4,2]),
  "R_4_perc"= c(Fisher_Grado_de_invasion_escala_Knosp_all$prop.col[5,2],
                Fisher_Grado_de_invasion_escala_Knosp_850k$prop.col[5,2],
                Fisher_Grado_de_invasion_escala_Knosp_935k$prop.col[5,2]),
  "p_value_true odds ratio is not equal to 1" = c(Fisher_Grado_de_invasion_escala_Knosp_all$fisher.ts$p.value,
                                          Fisher_Grado_de_invasion_escala_Knosp_850k$fisher.ts$p.value,
                                          Fisher_Grado_de_invasion_escala_Knosp_935k$fisher.ts$p.value))#, 
# "p_value_true odds ratio is less than 1" = c(Fisher_Grado_de_invasion_escala_Knosp_all$fisher.tl$p.value,
#                                       Fisher_Grado_de_invasion_escala_Knosp_850k$fisher.tl$p.value,
#                                       Fisher_Grado_de_invasion_escala_Knosp_935k$fisher.tl$p.value),
#  "p_value_true odds ratio is greater than 1" = c(Fisher_Grado_de_invasion_escala_Knosp_all$fisher.gt$p.value,
#                                          Fisher_Grado_de_invasion_escala_Knosp_850k$fisher.gt$p.value,
#                                          Fisher_Grado_de_invasion_escala_Knosp_935k$fisher.gt$p.value))
rownames(fisher_Grado_de_invasion_escala_Knosp_df) <- c("Grado_de_invasion_escala_Knosp_ALL",
                                                      "Grado_de_invasion_escala_Knosp_850k",
                                                      "Grado_de_invasion_escala_Knosp_935k")

write.table(fisher_Grado_de_invasion_escala_Knosp_df, file= "20240617_data_fisher_Grado_de_invasion_escala_Knosp.txt", 
            sep = "\t", quote = F, col.names = T, row.names = F)


#######################################################
## EDAD AL DIAGNOSTICO
library(tidyverse)
library(ggpubr)
library(rstatix)

# Check the Normal distribution
# Shapiro: From the output, the p-value is greater than the significance level 0.05 indicating that the distribution of the data are not significantly different from the normal distribution. In other words, we can assume the normality.
R_Norm_all <- shapiro_test(phenoData_all[which(phenoData_all$Group == "R"), "Edad_al_diagnostico"])
R_Norm_850k <- shapiro_test(phenoData_850k[which(phenoData_850k$Group == "R"), "Edad_al_diagnostico"])
R_Norm_935k <- shapiro_test(phenoData_935k[which(phenoData_935k$Group == "R"), "Edad_al_diagnostico"])
No_Norm_all <- shapiro_test(phenoData_all[which(phenoData_all$Group == "No"), "Edad_al_diagnostico"])
No_Norm_850k <- shapiro_test(phenoData_850k[which(phenoData_850k$Group == "No"), "Edad_al_diagnostico"])
No_Norm_935k <- shapiro_test(phenoData_935k[which(phenoData_935k$Group == "No"), "Edad_al_diagnostico"])

# Check the equality of variances
#This can be done using the Levene’s test. If the variances of groups are equal, the p-value should be greater than 0.05.
Var_all <- phenoData_all %>% levene_test(Edad_al_diagnostico ~ Group)
Var_850k <- phenoData_850k %>% levene_test(Edad_al_diagnostico ~ Group)
Var_935k <- phenoData_935k %>% levene_test(Edad_al_diagnostico ~ Group)
#If The p-value of the Levene’s test is significant, suggesting that there is a significant difference between the variances of the two groups. Therefore, we’ll use the Weltch t-test, which doesn’t assume the equality of the two variances.

# t-test
stat.test_all <- phenoData_all %>% t_test(Edad_al_diagnostico ~ Group) # %>% add_significance()
stat.test_850k <- phenoData_850k %>% t_test(Edad_al_diagnostico ~ Group) # %>% add_significance()
stat.test_935k <- phenoData_935k %>% t_test(Edad_al_diagnostico ~ Group) # %>% add_significance()

#Cohen’s d for Welch t-test
#The Welch test is a variant of t-test used when the equality of variance can’t be assumed. The effect size can be computed by dividing the mean difference between the groups by the “averaged” standard deviation.
phenoData_all %>% cohens_d(Edad_al_diagnostico ~ Group, var.equal = FALSE)


edad_means_sd_df <- data.frame(R_mean = c(mean(phenoData_all$Edad_al_diagnostico[which(phenoData_all$Group == "R")]),
                                     mean(phenoData_850k$Edad_al_diagnostico[which(phenoData_850k$Group == "R")]),
                                     mean(phenoData_935k$Edad_al_diagnostico[which(phenoData_935k$Group == "R")])),
                            R_sd = c(sd(phenoData_all$Edad_al_diagnostico[which(phenoData_all$Group == "R")]),
                                   sd(phenoData_850k$Edad_al_diagnostico[which(phenoData_850k$Group == "R")]),
                                   sd(phenoData_935k$Edad_al_diagnostico[which(phenoData_935k$Group == "R")])),
                            No_mean = c(mean(phenoData_all$Edad_al_diagnostico[which(phenoData_all$Group == "No")]),
                                       mean(phenoData_850k$Edad_al_diagnostico[which(phenoData_850k$Group == "No")]),
                                       mean(phenoData_935k$Edad_al_diagnostico[which(phenoData_935k$Group == "No")])),
                            No_sd = c(sd(phenoData_all$Edad_al_diagnostico[which(phenoData_all$Group == "No")]),
                                     sd(phenoData_850k$Edad_al_diagnostico[which(phenoData_850k$Group == "No")]),
                                     sd(phenoData_935k$Edad_al_diagnostico[which(phenoData_935k$Group == "No")])),
                            R_Shapiro_pVal = c(R_Norm_all$p.value, R_Norm_850k$p.value, R_Norm_935k$p.value),
                            No_Shapiro_pVal = c(No_Norm_all$p.value, No_Norm_850k$p.value, No_Norm_935k$p.value),
                            levene_pval = c(Var_all$p, Var_850k$p, Var_935k$p),
                            t_test_pVal = c(stat.test_all$p, stat.test_850k$p, stat.test_935k$p)) 

rownames(edad_means_sd_df) <- c("Edad_ALL","Edad_850k","Edad_935k")                   


write.table(edad_means_sd_df, file= "20240617_data_test_Edad_de_diagnostico.txt", 
            sep = "\t", quote = F, col.names = T, row.names = F)


#######################################################
## TAMAÑO DEL TUMOR MM

# Check the Normal distribution
# Shapiro: From the output, the p-value is greater than the significance level 0.05 indicating that the distribution of the data are not significantly different from the normal distribution. In other words, we can assume the normality.
R_Norm_all <- shapiro_test(phenoData_all[which(phenoData_all$Group == "R"), "Tamano_del_tumor_mm"])
R_Norm_850k <- shapiro_test(phenoData_850k[which(phenoData_850k$Group == "R"), "Tamano_del_tumor_mm"])
R_Norm_935k <- shapiro_test(phenoData_935k[which(phenoData_935k$Group == "R"), "Tamano_del_tumor_mm"])
No_Norm_all <- shapiro_test(phenoData_all[which(phenoData_all$Group == "No"), "Tamano_del_tumor_mm"])
No_Norm_850k <- shapiro_test(phenoData_850k[which(phenoData_850k$Group == "No"), "Tamano_del_tumor_mm"])
No_Norm_935k <- shapiro_test(phenoData_935k[which(phenoData_935k$Group == "No"), "Tamano_del_tumor_mm"])

# Check the equality of variances
#This can be done using the Levene’s test. If the variances of groups are equal, the p-value should be greater than 0.05.
Var_all <- phenoData_all %>% levene_test(Tamano_del_tumor_mm ~ Group)
Var_850k <- phenoData_850k %>% levene_test(Tamano_del_tumor_mm ~ Group)
Var_935k <- phenoData_935k %>% levene_test(Tamano_del_tumor_mm ~ Group)
#If The p-value of the Levene’s test is significant, suggesting that there is a significant difference between the variances of the two groups. Therefore, we’ll use the Weltch t-test, which doesn’t assume the equality of the two variances.

# t-texst
stat.test_all <- phenoData_all %>% t_test(Tamano_del_tumor_mm ~ Group) # %>% add_significance()
stat.test_850k <- phenoData_850k %>% t_test(Tamano_del_tumor_mm ~ Group) # %>% add_significance()
stat.test_935k <- phenoData_935k %>% t_test(Tamano_del_tumor_mm ~ Group) # %>% add_significance()
 

tamano_means_sd_df <- data.frame(R_mean = c(mean(phenoData_all$Tamano_del_tumor_mm[which(phenoData_all$Group == "R")], na.rm = TRUE),
                                          mean(phenoData_850k$Tamano_del_tumor_mm[which(phenoData_850k$Group == "R")], na.rm = TRUE),
                                          mean(phenoData_935k$Tamano_del_tumor_mm[which(phenoData_935k$Group == "R")], na.rm = TRUE)),
                               R_sd = c(sd(phenoData_all$Tamano_del_tumor_mm[which(phenoData_all$Group == "R")], na.rm = TRUE),
                                        sd(phenoData_850k$Tamano_del_tumor_mm[which(phenoData_850k$Group == "R")], na.rm = TRUE),
                                        sd(phenoData_935k$Tamano_del_tumor_mm[which(phenoData_935k$Group == "R")], na.rm = TRUE)),
                               No_mean = c(mean(phenoData_all$Tamano_del_tumor_mm[which(phenoData_all$Group == "No")], na.rm = TRUE),
                                           mean(phenoData_850k$Tamano_del_tumor_mm[which(phenoData_850k$Group == "No")], na.rm = TRUE),
                                           mean(phenoData_935k$Tamano_del_tumor_mm[which(phenoData_935k$Group == "No")], na.rm = TRUE)),
                               No_sd = c(sd(phenoData_all$Tamano_del_tumor_mm[which(phenoData_all$Group == "No")], na.rm = TRUE),
                                         sd(phenoData_850k$Tamano_del_tumor_mm[which(phenoData_850k$Group == "No")], na.rm = TRUE),
                                         sd(phenoData_935k$Tamano_del_tumor_mm[which(phenoData_935k$Group == "No")], na.rm = TRUE)),
                               R_Shapiro_pVal = c(R_Norm_all$p.value, R_Norm_850k$p.value, R_Norm_935k$p.value),
                               No_Shapiro_pVal = c(No_Norm_all$p.value, No_Norm_850k$p.value, No_Norm_935k$p.value),
                               levene_pval = c(Var_all$p, Var_850k$p, Var_935k$p),
                               t_test_pVal = c(stat.test_all$p, stat.test_850k$p, stat.test_935k$p)) 

rownames(tamano_means_sd_df) <- c("Tamano_ALL","Tamano_850k","Tamano_935k")                   
                                             
#Cohen’s d for Welch t-test
#The Welch test is a variant of t-test used when the equality of variance can’t be assumed. The effect size can be computed by dividing the mean difference between the groups by the “averaged” standard deviation.
phenoData_all %>% cohens_d(Tamano_del_tumor_mm ~ Group, var.equal = FALSE)


write.table(tamano_means_sd_df, file= "20240617_data_test_Tamano_del_tumor_mm.txt", 
            sep = "\t", quote = F, col.names = T, row.names = F)


#######################################################
## Ki67_porc_Num

# Check the Normal distribution
# Shapiro: From the output, the p-value is greater than the significance level 0.05 indicating that the distribution of the data are not significantly different from the normal distribution. In other words, we can assume the normality.
R_Norm_all <- shapiro_test(phenoData_all[which(phenoData_all$Group == "R"), "Ki67_porc_Num"])
R_Norm_850k <- shapiro_test(phenoData_850k[which(phenoData_850k$Group == "R"), "Ki67_porc_Num"])
R_Norm_935k <- shapiro_test(phenoData_935k[which(phenoData_935k$Group == "R"), "Ki67_porc_Num"])
No_Norm_all <- shapiro_test(phenoData_all[which(phenoData_all$Group == "No"), "Ki67_porc_Num"])
No_Norm_850k <- shapiro_test(phenoData_850k[which(phenoData_850k$Group == "No"), "Ki67_porc_Num"])
No_Norm_935k <- shapiro_test(phenoData_935k[which(phenoData_935k$Group == "No"), "Ki67_porc_Num"])

# Check the equality of variances
#This can be done using the Levene’s test. If the variances of groups are equal, the p-value should be greater than 0.05.
Var_all <- phenoData_all %>% levene_test(Ki67_porc_Num ~ Group)
Var_850k <- phenoData_850k %>% levene_test(Ki67_porc_Num ~ Group)
Var_935k <- phenoData_935k %>% levene_test(Ki67_porc_Num ~ Group)
#If The p-value of the Levene’s test is significant, suggesting that there is a significant difference between the variances of the two groups. Therefore, we’ll use the Weltch t-test, which doesn’t assume the equality of the two variances.

# t-texst
stat.test_all <- phenoData_all %>% t_test(Ki67_porc_Num ~ Group) # %>% add_significance()
stat.test_850k <- phenoData_850k %>% t_test(Ki67_porc_Num ~ Group) # %>% add_significance()
stat.test_935k <- phenoData_935k %>% t_test(Ki67_porc_Num ~ Group) # %>% add_significance()


Ki67_means_sd_df <- data.frame(R_mean = c(mean(phenoData_all$Ki67_porc_Num[which(phenoData_all$Group == "R")], na.rm = TRUE),
                                            mean(phenoData_850k$Ki67_porc_Num[which(phenoData_850k$Group == "R")], na.rm = TRUE),
                                            mean(phenoData_935k$Ki67_porc_Num[which(phenoData_935k$Group == "R")], na.rm = TRUE)),
                                 R_sd = c(sd(phenoData_all$Ki67_porc_Num[which(phenoData_all$Group == "R")], na.rm = TRUE),
                                          sd(phenoData_850k$Ki67_porc_Num[which(phenoData_850k$Group == "R")], na.rm = TRUE),
                                          sd(phenoData_935k$Ki67_porc_Num[which(phenoData_935k$Group == "R")], na.rm = TRUE)),
                                 No_mean = c(mean(phenoData_all$Ki67_porc_Num[which(phenoData_all$Group == "No")], na.rm = TRUE),
                                             mean(phenoData_850k$Ki67_porc_Num[which(phenoData_850k$Group == "No")], na.rm = TRUE),
                                             mean(phenoData_935k$Ki67_porc_Num[which(phenoData_935k$Group == "No")], na.rm = TRUE)),
                                 No_sd = c(sd(phenoData_all$Ki67_porc_Num[which(phenoData_all$Group == "No")], na.rm = TRUE),
                                           sd(phenoData_850k$Ki67_porc_Num[which(phenoData_850k$Group == "No")], na.rm = TRUE),
                                           sd(phenoData_935k$Ki67_porc_Num[which(phenoData_935k$Group == "No")], na.rm = TRUE)),
                                 R_Shapiro_pVal = c(R_Norm_all$p.value, R_Norm_850k$p.value, R_Norm_935k$p.value),
                                 No_Shapiro_pVal = c(No_Norm_all$p.value, No_Norm_850k$p.value, No_Norm_935k$p.value),
                                 levene_pval = c(Var_all$p, Var_850k$p, Var_935k$p)) 

rownames(Ki67_means_sd_df) <- c("Tamano_ALL","Tamano_850k","Tamano_935k")                   

#Cohen’s d for Welch t-test
#The Welch test is a variant of t-test used when the equality of variance can’t be assumed. The effect size can be computed by dividing the mean difference between the groups by the “averaged” standard deviation.
phenoData_all %>% cohens_d(Ki67_porc_Num ~ Group, var.equal = FALSE)

write.table(Ki67_means_sd_df, file= "20240617_data_test_Ki67_porc_Num.txt", 
            sep = "\t", quote = F, col.names = T, row.names = F)

