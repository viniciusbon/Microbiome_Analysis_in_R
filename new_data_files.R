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

# Importação dos dados
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

# Fase 1 - Pré-processamento e criação dos dataframes
df_ko <- data_storage$KO
df_AMR <- data_storage$AMR
df_VF <- data_storage$VF  # Corrigi de FV para VF
df_CAZy <- data_storage$CAZy
df_EC <- data_storage$EC
df_COGs <- data_storage$COGs
df_meta <- data_storage$meta

# Verificação dos formatos originais
cat("=== VERIFICAÇÃO DOS FORMATOS ORIGINAIS ===\n")
datasets <- list("KO" = df_ko, "AMR" = df_AMR, "VF" = df_VF, 
                 "CAZy" = df_CAZy, "EC" = df_EC, "COGs" = df_COGs)

for(name in names(datasets)) {
  if(!is.null(datasets[[name]])) {
    dims <- dim(datasets[[name]])
    cat(paste(name, ":", dims[1], "linhas x", dims[2], "colunas\n"))
    cat("Primeiras colunas:", paste(names(datasets[[name]])[1:min(3, ncol(datasets[[name]]))], collapse = ", "), "\n\n")
  }
}

# FASE 2 - TRANSPOSIÇÃO E PADRONIZAÇÃO

# Função para transpor e padronizar
transpose_functional_data <- function(data, feature_id_name) {
  # Salvar os IDs das features (primeira coluna)
  feature_ids <- data[[1]]
  
  # Salvar os nomes das amostras (nomes das colunas, exceto a primeira)
  sample_names <- names(data)[-1]
  
  # Extrair apenas os dados numéricos (excluindo primeira coluna)
  count_matrix <- data[, -1]
  
  # Transpor a matriz
  transposed_matrix <- t(count_matrix)
  
  # Converter para dataframe
  transposed_df <- as.data.frame(transposed_matrix)
  
  # Definir nomes das colunas como os feature IDs
  names(transposed_df) <- feature_ids
  
  # Adicionar coluna Sample_ID
  transposed_df <- transposed_df %>%
    mutate(Sample_ID = sample_names, .before = 1)
  
  return(transposed_df)
}

# TRANSPOSIÇÃO DOS DATASETS ESPECIFICADOS

cat("=== INICIANDO TRANSPOSIÇÕES ===\n\n")

# 1. Transpor df_ko
cat("1. Transpondendo df_ko...\n")
if(!is.null(df_ko)) {
  df_ko_transposed <- transpose_functional_data(df_ko, "ko")
  cat("✓ df_ko transposto:", nrow(df_ko_transposed), "linhas x", ncol(df_ko_transposed), "colunas\n")
  df_ko <- df_ko_transposed
} else {
  cat("⚠ df_ko não encontrado\n")
}

# 2. Transpor df_CAZy
cat("2. Transpondendo df_CAZy...\n")
if(!is.null(df_CAZy)) {
  df_CAZy_transposed <- transpose_functional_data(df_CAZy, "cazy")
  cat("✓ df_CAZy transposto:", nrow(df_CAZy_transposed), "linhas x", ncol(df_CAZy_transposed), "colunas\n")
  df_CAZy <- df_CAZy_transposed
} else {
  cat("⚠ df_CAZy não encontrado\n")
}

# 3. Transpor df_COGs
cat("3. Transpondendo df_COGs...\n")
if(!is.null(df_COGs)) {
  df_COGs_transposed <- transpose_functional_data(df_COGs, "cog")
  cat("✓ df_COGs transposto:", nrow(df_COGs_transposed), "linhas x", ncol(df_COGs_transposed), "colunas\n")
  df_COGs <- df_COGs_transposed
} else {
  cat("⚠ df_COGs não encontrado\n")
}

# 4. Transpor df_EC
cat("4. Transpondendo df_EC...\n")
if(!is.null(df_EC)) {
  df_EC_transposed <- transpose_functional_data(df_EC, "ec")
  cat("✓ df_EC transposto:", nrow(df_EC_transposed), "linhas x", ncol(df_EC_transposed), "colunas\n")
  df_EC <- df_EC_transposed
} else {
  cat("⚠ df_EC não encontrado\n")
}

# VERIFICAÇÃO FINAL DOS FORMATOS
cat("\n=== VERIFICAÇÃO FINAL DOS FORMATOS PADRONIZADOS ===\n")
datasets_final <- list("KO" = df_ko, "AMR" = df_AMR, "VF" = df_VF, 
                       "CAZy" = df_CAZy, "EC" = df_EC, "COGs" = df_COGs)

for(name in names(datasets_final)) {
  if(!is.null(datasets_final[[name]])) {
    dims <- dim(datasets_final[[name]])
    first_col <- names(datasets_final[[name]])[1]
    cat(paste(name, ":", dims[1], "linhas x", dims[2], "colunas"))
    cat(paste(" | Primeira coluna:", first_col, "\n"))
  }
}

# Preview dos dados transpostos
cat("\n=== PREVIEW DOS DADOS PADRONIZADOS ===\n")
for(name in names(datasets_final)) {
  if(!is.null(datasets_final[[name]])) {
    cat(paste("\n", name, "- Primeiras 3 linhas e 5 colunas:\n"))
    print(datasets_final[[name]][1:min(3, nrow(datasets_final[[name]])), 
                                 1:min(5, ncol(datasets_final[[name]]))])
  }
}

cat("\n✅ TRANSPOSIÇÃO CONCLUÍDA! Todos os datasets agora estão padronizados:\n")
cat("- Amostras como linhas\n")
cat("- Features funcionais como colunas\n") 
cat("- Primeira coluna = Sample_ID\n")
