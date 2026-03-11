#Pipeline to analyse: 
#fun_VF_counts
#fun_KO_counts
#fun_EC_count 
#fun_COGs_count
#fun_CAZy_counts, 
#fun_AMR_counts.tsv

#import libs
# Core preprocessing packages
library(phyloseq)        # Framework microbioma
library(DESeq2)          # Normalização e análise diferencial
library(edgeR)           # TMM normalization
library(compositions)    # CLR transformation
library(zCompositions)   # Tratamento de zeros
library(dplyr)           # Manipulação de dados
library(tibble)          # Data frames modernos
library(readr)           # Leitura de arquivos
library(stringr)         # Manipulação de strings
library(ggplot2)         # Visualizações básicas
library(corrplot)        # Gráficos de correlação
library(VennDiagram)     # Diagramas de Venn
library(readr)
library(readxl)

#Fase 0 - Importing the datasets

#Set the data dir where your samples are located

data_dir <- "C:/Users/dti-/Documents/ANH/Simphyome_Analysis/Hav1801" 

#Files dictionary
patterns <- c(
  "AMR"   = "AMR_counts",
  "CAZy"  = "CAZy_counts",
  "COGs"  = "COGs_count",
  "EC"    = "EC_count",
  "KO"    = "KO_counts",
  "VF"    = "VF_counts",
  "meta"  = "meta"
)

#data storage object
data_storage <- list()

# Search and Import
for (key in names(patterns)) {
  # Search for all the files using the dictionary
  file_found <- list.files(path = data_dir, pattern = patterns[key], full.names = TRUE)
  
  if (length(file_found) == 0) {
    warning(paste("No files found for:", key))
    next
  } else if (length(file_found) > 1) {
    message(paste("Files found for:", key, "using:", basename(file_found[1])))
    file_found <- file_found[1]
  }
  
  if (grepl("\\.xlsx$", file_found)) {
    data_storage[[key]] <- read_excel(file_found)
  } else {
    data_storage[[key]] <- read_delim(file_found, delim = "\t", show_col_types = FALSE)
  }
  
  message(paste("Success = ", key))
}

