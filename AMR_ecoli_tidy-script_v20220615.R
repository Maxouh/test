# On télécharge les packages nécessaires
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# On définit l'environnement de travail
rm(list = ls())
setwd("/Users/maxou/Desktop/PhD_ecoli_tmp")


# Formattage des données
# Etape 1 : On télécharge le fichier (matrice de presence/absence fourni par abricate --summary)
m <- read.csv("resfinder_matrix", sep = "\t", header = TRUE)

# Etape 2 : On tidy les données en 3 colonnes
# colonnes : "Assembly.barcode","bla_gene", "presence_absence"
colnames(m)[1] <- "Assembly.barcode"

# On coupe la partie des noms de souches contenant le chemin d'accès utilisé par abricate
m$Assembly.barcode <- sub("contigs_included/", "", m$Assembly.barcode)

  # On sélectionne que les gènes bla
  m <- m[,c(1,grep("bla", colnames(m)))]
  # On transforme les valeurs "." en absence
  m[m == "."] <- "absence"
  # On créé un tableau n à partir de m en omettant la colonne "Assembly.barcode"
  # On dit que tout ce qui est différente de "absence" -> "presence"
  # D'où le fait de retirer le nom des souches (sinon elles sont toutes appelées "presence")
  n <- m[,-1]
  n[n != "absence"] <- "presence"
  
  # On simplifie le nom des blaGENE en retirant le "_1" via la fonction X <- gsub("\\_1","",X)
  colnames(m)[-1] <- gsub("\\_1","",colnames(m)[-1])

  # On colle la colonne Assembly.barcode avec la nouvelle matrice de presence_absence et on supprime n
  m <- cbind(m$Assembly.barcode,n)
  rm(n)
  
  # On renomme la colonne "m$Assembly.barcode" en juste "Assembly.barcode"
  colnames(m)[1] <- "Assembly.barcode"

  # On simplifie le nom des blaGENE en retirant le "_1" via la fonction X <- gsub("\\_*","",X)
  colnames(m)[-1] <- gsub("\\_1","",colnames(m)[-1])
  
  # On veut améliorer l'écriture des bla_gene ==> blaCTX.M.3_1 ->  blaCTX.M.3 -> blaCTX.M_3 -> blaCTX-M_3
  colnames(m)[-1] <- gsub("\\_1","",colnames(m)[-1])
  colnames(m)[-1] <- gsub("blaCTX.M","blaCTX-M",colnames(m)[-1])
  colnames(m)[-1] <- gsub("[.]","_", colnames(m)[-1])
  
  # On pivot_longer toute la matrice pour n'avoir que 3 colonnes "Assembly;barcode", "bla_gene" et "presence_absence"
  m <- tibble(m)
  m <- m %>% pivot_longer(!Assembly.barcode, names_to = "bla_gene", values_to = "presence_absence")
  m <- m %>% filter(presence_absence == "presence")   # On ne prend que les genes effectivement present dans chacune que des souches
  # m <- m %>% filter(bla_gene != "blaZ_79")        # pas de blaZ_79 dans les souches... (abricate --summary a retrouvé 997 hits de blaZ_79 -> connerie +++)
  

# Etape 3 : on dit qui est BLSE et qui ne l'est 
  # Etape 1 : on écrit le fichier bldb_esbl temporaire (avec uniquement les bla_genes retrouvés dans l'étude)
  write.csv(sort(unique(m$bla_gene)), "bldb_esbl.csv", row.names = FALSE, quote = FALSE)
  
   "" [...] on créé le fichier bla_db à partir de bldb_esbl.csv et de bldb pour donner les informations ESBL/EPC/bla_gene_group (3 nouvelles colonnes) ""
                                                                                                                                
  # Etape 2 : On rempli le fichier csv sur excel MANUELLEMENT et on l'enregistre au nom "bla_db.csv"
  
  # Etape 3 : On récupère le fichier avec les ajouts
  d <- read.csv("bla_db.csv", header = TRUE, sep = ";")
  d <- tibble(d)
  
  # On fusionne les deux tableaux (m et d) en prenant la variable bla_gene pour base
  df <- left_join(m,d, by = "bla_gene")
  write.table(df, "amr_metadata.csv", row.names = FALSE, quote = FALSE)
  

  getwd()
  list.files()
  getwd()





