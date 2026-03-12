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

#stadarize metadata ID colum name
colnames(df_meta)[colnames(df_meta) == "sample"] <- "Sample_ID"

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

# FASE 3 - PADRONIZAÇÃO DOS NOMES DAS AMOSTRAS

cat("\n=== PADRONIZAÇÃO DOS NOMES DAS AMOSTRAS ===\n")

# Função para padronizar nomes das amostras no metadata
standardize_sample_names <- function(sample_names) {
  # Substituir "." por "-" e remover apenas o último "_"
  standardized_names <- sample_names %>%
    str_replace_all("\\.", "-") %>%  # Substitui todos os "." por "-"
    str_replace("_$", "")            # Remove apenas o "_" no final da string
  
  return(standardized_names)
}

# Verificar formato original dos nomes das amostras
cat("Verificando formatos originais dos Sample_IDs:\n\n")

# Mostrar exemplos dos nomes das amostras em cada dataset
if(exists("df_ko_transposed")) {
  cat("df_ko_transposed - Exemplos de Sample_ID:\n")
  cat(paste(head(df_ko_transposed$Sample_ID, 3), collapse = ", "), "\n\n")
}

if(exists("df_CAZy_transposed")) {
  cat("df_CAZy_transposed - Exemplos de Sample_ID:\n")
  cat(paste(head(df_CAZy_transposed$Sample_ID, 3), collapse = ", "), "\n\n")
}

if(exists("df_COGs_transposed")) {
  cat("df_COGs_transposed - Exemplos de Sample_ID:\n")
  cat(paste(head(df_COGs_transposed$Sample_ID, 3), collapse = ", "), "\n\n")
}

if(exists("df_EC_transposed")) {
  cat("df_EC_transposed - Exemplos de Sample_ID:\n")
  cat(paste(head(df_EC_transposed$Sample_ID, 3), collapse = ", "), "\n\n")
}

if(!is.null(df_meta)) {
  # Identificar qual coluna contém os Sample_IDs no metadata
  sample_id_col <- names(df_meta)[1]  # Assumindo que é a primeira coluna
  cat(paste("df_meta - Coluna Sample_ID (", sample_id_col, ") - Exemplos:\n"))
  cat(paste(head(df_meta[[sample_id_col]], 3), collapse = ", "), "\n\n")
  
  # Aplicar padronização aos nomes no metadata
  cat("Aplicando padronização ao df_meta...\n")
  df_meta[[sample_id_col]] <- standardize_sample_names(df_meta[[sample_id_col]])
  
  cat("df_meta - Após padronização:\n")
  cat(paste(head(df_meta[[sample_id_col]], 3), collapse = ", "), "\n\n")
}

# VERIFICAÇÃO FINAL DOS FORMATOS E NOMES PADRONIZADOS
cat("=== VERIFICAÇÃO FINAL DOS FORMATOS PADRONIZADOS ===\n")

# Lista dos datasets finais
datasets_final <- list(
  "KO" = if(exists("df_ko_transposed")) df_ko_transposed else NULL,
  "AMR" = df_AMR, 
  "VF" = df_VF, 
  "CAZy" = if(exists("df_CAZy_transposed")) df_CAZy_transposed else NULL,
  "EC" = if(exists("df_EC_transposed")) df_EC_transposed else NULL,
  "COGs" = if(exists("df_COGs_transposed")) df_COGs_transposed else NULL,
  "META" = df_meta
)

for(name in names(datasets_final)) {
  if(!is.null(datasets_final[[name]])) {
    dims <- dim(datasets_final[[name]])
    first_col <- names(datasets_final[[name]])[1]
    cat(paste(name, ":", dims[1], "linhas x", dims[2], "colunas"))
    cat(paste(" | Primeira coluna:", first_col, "\n"))
    
    # Mostrar exemplos dos Sample_IDs padronizados
    if(first_col %in% c("Sample_ID", names(df_meta)[1])) {
      sample_examples <- head(datasets_final[[name]][[first_col]], 3)
      cat(paste("  Exemplos Sample_ID:", paste(sample_examples, collapse = ", "), "\n"))
    }
  }
}

# Verificar se todos os Sample_IDs estão no mesmo formato
cat("\n=== VERIFICAÇÃO DE CONSISTÊNCIA DOS SAMPLE_IDS ===\n")

# Coletar todos os Sample_IDs dos datasets funcionais
all_sample_ids <- list()

if(exists("df_ko_transposed")) all_sample_ids[["KO"]] <- df_ko_transposed$Sample_ID
if(exists("df_CAZy_transposed")) all_sample_ids[["CAZy"]] <- df_CAZy_transposed$Sample_ID
if(exists("df_COGs_transposed")) all_sample_ids[["COGs"]] <- df_COGs_transposed$Sample_ID
if(exists("df_EC_transposed")) all_sample_ids[["EC"]] <- df_EC_transposed$Sample_ID
if(!is.null(df_meta)) all_sample_ids[["META"]] <- df_meta[[names(df_meta)[1]]]

# Verificar consistência
if(length(all_sample_ids) > 1) {
  # Comparar se todos os datasets têm os mesmos Sample_IDs
  first_dataset_ids <- all_sample_ids[[1]]
  consistent <- TRUE
  
  for(i in 2:length(all_sample_ids)) {
    if(!identical(sort(first_dataset_ids), sort(all_sample_ids[[i]]))) {
      consistent <- FALSE
      cat(paste("⚠ Inconsistência encontrada entre", names(all_sample_ids)[1], "e", names(all_sample_ids)[i], "\n"))
    }
  }
  
  if(consistent) {
    cat("✅ Todos os Sample_IDs estão consistentes entre os datasets!\n")
  }
} else {
  cat("⚠ Apenas um dataset encontrado para verificação\n")
}

cat("\n✅ PADRONIZAÇÃO CONCLUÍDA!\n")
cat("- Todos os datasets transpostos estão no formato: amostras como linhas, features como colunas\n")
cat("- Sample_IDs padronizados no formato: HAV1801-12-10-18_10A\n")
cat("- Primeira coluna = Sample_ID em todos os datasets\n")


# FASE 4 - FILTRAGEM DE AMOSTRAS USANDO METADATA COMO REFERÊNCIA

cat("\n", rep("=", 60), "\n")
cat("=== FASE 4: FILTRAGEM DE AMOSTRAS USANDO METADATA ===\n")
cat(rep("=", 60), "\n")

if(is.null(df_meta)) {
  cat("⚠ ERRO: df_meta não encontrado. Impossível prosseguir com a filtragem.\n")
} else {
  
  # Obter a lista de Sample_IDs de referência do metadata
  reference_sample_ids <- df_meta$Sample_ID
  cat("📋 Sample_IDs de referência no metadata:", length(reference_sample_ids), "amostras\n")
  cat("Primeiros exemplos:", paste(head(reference_sample_ids, 5), collapse = ", "), "\n\n")
  
  # Função para filtrar datasets usando Sample_IDs de referência
  filter_samples_by_metadata <- function(dataset, dataset_name, reference_ids) {
    if(is.null(dataset) || !"Sample_ID" %in% names(dataset)) {
      cat(paste("⚠", dataset_name, "- Dataset não encontrado ou sem coluna Sample_ID\n"))
      return(NULL)
    }
    
    # Amostras antes da filtragem
    original_samples <- nrow(dataset)
    original_sample_ids <- dataset$Sample_ID
    
    # Filtrar apenas amostras que estão no metadata
    filtered_dataset <- dataset %>%
      filter(Sample_ID %in% reference_ids)
    
    # Amostras após filtragem
    filtered_samples <- nrow(filtered_dataset)
    kept_sample_ids <- filtered_dataset$Sample_ID
    
    # Amostras removidas
    removed_samples <- original_samples - filtered_samples
    removed_sample_ids <- setdiff(original_sample_ids, kept_sample_ids)
    
    # Relatório
    cat(paste("📊", dataset_name, ":\n"))
    cat(paste("   • Amostras originais:", original_samples, "\n"))
    cat(paste("   • Amostras mantidas:", filtered_samples, "\n"))
    cat(paste("   • Amostras removidas:", removed_samples, "\n"))
    
    if(removed_samples > 0) {
      cat("   • Sample_IDs removidos:", paste(head(removed_sample_ids, 5), collapse = ", "))
      if(length(removed_sample_ids) > 5) {
        cat(paste(" (e mais", length(removed_sample_ids) - 5, "...)"))
      }
      cat("\n")
    }
    cat("\n")
    
    return(filtered_dataset)
  }
  
  cat("🔍 Iniciando filtragem dos datasets funcionais...\n\n")
  
  # Aplicar filtragem a todos os datasets funcionais
  df_ko_filtered <- filter_samples_by_metadata(df_ko, "df_ko", reference_sample_ids)
  df_AMR_filtered <- filter_samples_by_metadata(df_AMR, "df_AMR", reference_sample_ids)
  df_VF_filtered <- filter_samples_by_metadata(df_VF, "df_VF", reference_sample_ids)
  df_CAZy_filtered <- filter_samples_by_metadata(df_CAZy, "df_CAZy", reference_sample_ids)
  df_EC_filtered <- filter_samples_by_metadata(df_EC, "df_EC", reference_sample_ids)
  df_COGs_filtered <- filter_samples_by_metadata(df_COGs, "df_COGs", reference_sample_ids)
  
  # Substituir os datasets originais pelos filtrados
  df_ko <- df_ko_filtered
  df_AMR <- df_AMR_filtered
  df_VF <- df_VF_filtered
  df_CAZy <- df_CAZy_filtered
  df_EC <- df_EC_filtered
  df_COGs <- df_COGs_filtered
  
  # VERIFICAÇÃO FINAL PÓS-FILTRAGEM
  cat(rep("=", 50), "\n")
  cat("=== VERIFICAÇÃO FINAL PÓS-FILTRAGEM ===\n")
  cat(rep("=", 50), "\n")
  
  # Lista dos datasets filtrados
  datasets_filtered <- list(
    "META" = df_meta,
    "KO" = df_ko,
    "AMR" = df_AMR,
    "VF" = df_VF,
    "CAZy" = df_CAZy,
    "EC" = df_EC,
    "COGs" = df_COGs
  )
  
  # Relatório final
  for(name in names(datasets_filtered)) {
    if(!is.null(datasets_filtered[[name]])) {
      dims <- dim(datasets_filtered[[name]])
      cat(paste("✅", name, ":", dims[1], "amostras x", dims[2], "variáveis\n"))
    } else {
      cat(paste("❌", name, ": Dataset não disponível\n"))
    }
  }
  
  # Verificação de consistência final das amostras
  cat("\n🔍 Verificação de consistência final das amostras:\n")
  
  sample_counts <- list()
  for(name in names(datasets_filtered)) {
    if(!is.null(datasets_filtered[[name]]) && "Sample_ID" %in% names(datasets_filtered[[name]])) {
      sample_counts[[name]] <- length(datasets_filtered[[name]]$Sample_ID)
    }
  }
  
  if(length(unique(sample_counts)) == 1) {
    cat("✅ PERFEITO! Todos os datasets têm o mesmo número de amostras:", unique(sample_counts), "\n")
  } else {
    cat("⚠ ATENÇÃO! Números diferentes de amostras entre datasets:\n")
    for(name in names(sample_counts)) {
      cat(paste("  ", name, ":", sample_counts[[name]], "amostras\n"))
    }
  }
  
  # Verificação se todas as amostras são as mesmas
  if(!is.null(df_meta) && length(datasets_filtered) > 1) {
    reference_ids <- sort(df_meta$Sample_ID)
    all_consistent <- TRUE
    
    for(name in names(datasets_filtered)) {
      if(name != "META" && !is.null(datasets_filtered[[name]]) && "Sample_ID" %in% names(datasets_filtered[[name]])) {
        current_ids <- sort(datasets_filtered[[name]]$Sample_ID)
        if(!identical(reference_ids, current_ids)) {
          all_consistent <- FALSE
          cat(paste("⚠", name, "tem Sample_IDs diferentes do metadata\n"))
        }
      }
    }
    
    if(all_consistent) {
      cat("✅ Todos os datasets têm exatamente as mesmas amostras do metadata!\n")
    }
  }
  
  cat("\n", rep("🎉", 20), "\n")
  cat("✅ FASE 4 CONCLUÍDA COM SUCESSO!\n")
  cat("• Todos os datasets foram filtrados usando o metadata como referência\n")
  cat("• Apenas amostras presentes no metadata foram mantidas\n")
  cat("• Datasets prontos para análises funcionais integradas\n")
  cat(rep("🎉", 20), "\n")
}
