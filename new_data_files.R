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
selected_groups <- c("control", "glycodex")

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


library(dplyr)

# FASE 5 - JOIN DAS COLUNAS METADATA

cat("\n", rep("=", 60), "\n")
cat("=== FASE 5: ADICIONANDO COLUNAS AGE E TREATMENT ===\n")
cat(rep("=", 60), "\n")

# 1. Preparar subset do metadata apenas com as colunas necessárias
metadata_subset <- df_meta %>%
  select(Sample_ID, age, treatment)

cat("📋 Metadata subset preparado:\n")
print(head(metadata_subset, 3))
cat("\n")

# 2. Lista dos dataframes para fazer join
datasets_to_join <- list(
  "df_ko" = df_ko,
  "df_AMR" = df_AMR,
  "df_VF" = df_VF,
  "df_CAZy" = df_CAZy,
  "df_EC" = df_EC,
  "df_COGs" = df_COGs
)

# 3. Fazer left_join para cada dataset
cat("🔗 Fazendo join das colunas metadata:\n\n")

for(dataset_name in names(datasets_to_join)) {
  
  current_dataset <- datasets_to_join[[dataset_name]]
  
  if(!is.null(current_dataset)) {
    
    cat(paste("Processando", dataset_name, "...\n"))
    
    # Dimensões antes do join
    original_dims <- dim(current_dataset)
    
    # Fazer left_join
    dataset_with_metadata <- current_dataset %>%
      left_join(metadata_subset, by = "Sample_ID")
    
    # Dimensões após join
    new_dims <- dim(dataset_with_metadata)
    
    # Reordenar colunas: Sample_ID, age, treatment, depois resto
    cols_order <- c("Sample_ID", "age", "treatment", 
                    setdiff(names(dataset_with_metadata), c("Sample_ID", "age", "treatment")))
    
    dataset_with_metadata <- dataset_with_metadata[, cols_order]
    
    # Substituir o dataset original
    assign(dataset_name, dataset_with_metadata, envir = .GlobalEnv)
    
    # Relatório
    cat(paste("   ✅ Dimensões:", original_dims[1], "x", original_dims[2], "→", 
              new_dims[1], "x", new_dims[2], "\n"))
    cat(paste("   ✅ Colunas adicionadas:", new_dims[2] - original_dims[2], "\n"))
    
    # Verificar NAs
    na_age <- sum(is.na(dataset_with_metadata$age))
    na_treatment <- sum(is.na(dataset_with_metadata$treatment))
    cat(paste("   ✅ NAs em age:", na_age, "| NAs em treatment:", na_treatment, "\n"))
    
    cat("\n")
    
  } else {
    cat(paste("   ⚠", dataset_name, "não encontrado\n\n"))
  }
}

# 4. Verificação final
cat(rep("=", 50), "\n")
cat("=== VERIFICAÇÃO FINAL DOS JOINS ===\n")
cat(rep("=", 50), "\n")

final_datasets <- list(
  "df_meta" = df_meta,
  "df_ko" = df_ko,
  "df_AMR" = df_AMR,
  "df_VF" = df_VF,
  "df_CAZy" = df_CAZy,
  "df_EC" = df_EC,
  "df_COGs" = df_COGs
)

for(name in names(final_datasets)) {
  if(!is.null(final_datasets[[name]])) {
    dims <- dim(final_datasets[[name]])
    cols <- names(final_datasets[[name]])
    
    cat(paste("📊", name, ":", dims[1], "linhas x", dims[2], "colunas\n"))
    cat(paste("   Primeiras colunas:", paste(cols[1:min(6, length(cols))], collapse = ", "), "\n"))
    
    # Verificar se tem age e treatment
    has_age <- "age" %in% cols
    has_treatment <- "treatment" %in% cols
    cat(paste("   Age:", ifelse(has_age, "✅", "❌"), "| Treatment:", ifelse(has_treatment, "✅", "❌"), "\n\n"))
  }
}

# FASE 5 - JOIN SIMPLIFICADO DAS COLUNAS METADATA

cat("\n", rep("=", 60), "\n")
cat("=== FASE 5: ADICIONANDO AGE E TREATMENT (MÉTODO OTIMIZADO) ===\n")
cat(rep("=", 60), "\n")

library(dplyr)

# 1. Lista nomeada dos dataframes para facilitar identificação
dfs <- list(
  df_ko = df_ko,
  df_AMR = df_AMR,
  df_VF = df_VF,
  df_CAZy = df_CAZy,
  df_EC = df_EC,
  df_COGs = df_COGs
)

cat("📋 Processando", length(dfs), "datasets...\n")
cat("Metadata subset:", ncol(df_meta[, c("Sample_ID", "age", "treatment")]), "colunas\n\n")

# 2. Aplicar left_join com relatório
dfs_joined <- lapply(names(dfs), function(name) {
  
  df <- dfs[[name]]
  
  if(!is.null(df)) {
    
    cat(paste("🔗 Processando", name, "...\n"))
    
    # Dimensões antes
    original_dims <- dim(df)
    
    # Join
    df_with_metadata <- left_join(df, 
                                  df_meta[, c("Sample_ID", "age", "treatment")], 
                                  by = "Sample_ID")
    
    # Reordenar colunas: Sample_ID, age, treatment, depois resto
    cols_order <- c("Sample_ID", "age", "treatment", 
                    setdiff(names(df_with_metadata), c("Sample_ID", "age", "treatment")))
    df_with_metadata <- df_with_metadata[, cols_order]
    
    # Dimensões após
    new_dims <- dim(df_with_metadata)
    
    # Relatório
    cat(paste("   ✅", original_dims[1], "x", original_dims[2], "→", 
              new_dims[1], "x", new_dims[2], 
              "(+", new_dims[2] - original_dims[2], "colunas)\n"))
    
    # Verificar NAs
    na_age <- sum(is.na(df_with_metadata$age))
    na_treatment <- sum(is.na(df_with_metadata$treatment))
    cat(paste("   📊 NAs - Age:", na_age, "| Treatment:", na_treatment, "\n\n"))
    
    return(df_with_metadata)
    
  } else {
    cat(paste("   ⚠", name, "não encontrado\n\n"))
    return(NULL)
  }
})

# 3. Nomear a lista resultado
names(dfs_joined) <- names(dfs)

# 4. Substituir os dataframes originais no ambiente global
list2env(dfs_joined, envir = .GlobalEnv)

# 5. Verificação final compacta
cat(rep("=", 50), "\n")
cat("=== VERIFICAÇÃO FINAL ===\n")
cat(rep("=", 50), "\n")

verification_dfs <- list(df_ko, df_AMR, df_VF, df_CAZy, df_EC, df_COGs)
names(verification_dfs) <- c("df_ko", "df_AMR", "df_VF", "df_CAZy", "df_EC", "df_COGs")

for(name in names(verification_dfs)) {
  df <- verification_dfs[[name]]
  if(!is.null(df)) {
    dims <- dim(df)
    has_age <- "age" %in% names(df)
    has_treatment <- "treatment" %in% names(df)
    
    cat(paste("📊", name, ":", dims[1], "x", dims[2], 
              "| Age:", ifelse(has_age, "✅", "❌"),
              "| Treatment:", ifelse(has_treatment, "✅", "❌"), "\n"))
  }
}



# FASE 6 - FILTRO POR GRUPOS (CONTROL & GLYCODEX)

cat("\n", rep("=", 60), "\n")
cat("=== FASE 6: FILTRO POR GRUPOS (CONTROL & GLYCODEX) ===\n")
cat(rep("=", 60), "\n")

library(dplyr)

# 1. Grupos selecionados

cat("🎯 Grupos selecionados:", paste(selected_groups, collapse = ", "), "\n")

# 2. Verificar distribuição atual dos grupos
cat("\n📊 Distribuição atual dos grupos:\n")
treatment_counts <- table(df_meta$treatment)
print(treatment_counts)
cat("\n")

# 3. Lista dos dataframes para filtrar
dfs_to_filter <- list(
  df_ko = df_ko,
  df_AMR = df_AMR,
  df_VF = df_VF,
  df_CAZy = df_CAZy,
  df_EC = df_EC,
  df_COGs = df_COGs
)

# 4. Aplicar filtro com lapply (método eficiente)
cat("🔍 Aplicando filtro por grupos...\n\n")

dfs_grouped <- lapply(names(dfs_to_filter), function(name) {
  
  df <- dfs_to_filter[[name]]
  
  if(!is.null(df) && "treatment" %in% names(df)) {
    
    cat(paste("🔗 Filtrando", name, "...\n"))
    
    # Dimensões antes do filtro
    original_dims <- dim(df)
    original_groups <- table(df$treatment)
    
    # Filtrar por grupos selecionados
    df_filtered <- df %>%
      filter(treatment %in% selected_groups)
    
    # Dimensões após filtro
    new_dims <- dim(df_filtered)
    new_groups <- table(df_filtered$treatment)
    
    # Relatório
    cat(paste("   📏 Dimensões:", original_dims[1], "x", original_dims[2], "→", 
              new_dims[1], "x", new_dims[2], "\n"))
    cat(paste("   📊 Amostras removidas:", original_dims[1] - new_dims[1], "\n"))
    cat("   📈 Grupos finais:\n")
    print(new_groups)
    cat("\n")
    
    return(df_filtered)
    
  } else {
    cat(paste("   ⚠", name, "não encontrado ou sem coluna treatment\n\n"))
    return(NULL)
  }
})

# 5. Nomear e criar novos dataframes com sufixo "_grouped"
names(dfs_grouped) <- paste0(names(dfs_to_filter), "_grouped")

# 6. Criar as variáveis no ambiente global
list2env(dfs_grouped, envir = .GlobalEnv)

# 7. Verificação final
cat(rep("=", 50), "\n")
cat("=== VERIFICAÇÃO FINAL - DATASETS AGRUPADOS ===\n")
cat(rep("=", 50), "\n")

grouped_datasets <- list(
  df_ko_grouped, df_AMR_grouped, df_VF_grouped, 
  df_CAZy_grouped, df_EC_grouped, df_COGs_grouped
)
names(grouped_datasets) <- names(dfs_grouped)

for(name in names(grouped_datasets)) {
  df <- grouped_datasets[[name]]
  if(!is.null(df)) {
    dims <- dim(df)
    groups <- table(df$treatment)
    
    cat(paste("📊", name, ":", dims[1], "amostras x", dims[2], "variáveis\n"))
    cat(paste("   🎯 Grupos:", paste(names(groups), collapse = ", "), 
              "| Distribuição:", paste(groups, collapse = ", "), "\n\n"))
  }
}

cat("✅ FILTRO POR GRUPOS CONCLUÍDO!\n")
cat("• Novos datasets criados: df_*_grouped\n")
cat("• Apenas amostras 'control' e 'glycodex' mantidas\n")
cat("• Prontos para análises comparativas entre grupos!\n")



#Sample SWAP check
# PCA PARA VERIFICAÇÃO DE TROCA DE AMOSTRAS (GRUPOS SELECIONADOS)

cat("\n", rep("=", 60), "\n")
cat("=== PCA: VERIFICAÇÃO DE TROCA DE AMOSTRAS (CONTROL vs GLYCODEX) ===\n")
cat(rep("=", 60), "\n")

library(dplyr)
library(ggplot2)
library(ggrepel)
library(compositions)  # Para CLR transformation
library(tibble)

# Função para realizar PCA em dataset funcional
perform_functional_pca <- function(dataset, dataset_name, group_var = "treatment") {
  
  if(is.null(dataset) || !group_var %in% names(dataset)) {
    cat(paste("⚠", dataset_name, "não disponível ou sem coluna", group_var, "\n"))
    return(NULL)
  }
  
  cat(paste("🔍 Processando PCA para", dataset_name, "...\n"))
  
  # Verificar grupos disponíveis
  available_groups <- unique(dataset[[group_var]])
  cat(paste("   🎯 Grupos encontrados:", paste(available_groups, collapse = ", "), "\n"))
  
  # 1. Separar metadata das features funcionais
  metadata_cols <- c("Sample_ID", "age", "treatment")
  feature_cols <- setdiff(names(dataset), metadata_cols)
  
  # Extrair metadata e features
  metadata_subset <- dataset[, metadata_cols]
  features_data <- dataset[, feature_cols]
  
  # 2. Remover features com variância zero ou com muitos zeros
  # Filtrar colunas com pelo menos alguma variabilidade
  feature_vars <- apply(features_data, 2, var, na.rm = TRUE)
  valid_features <- feature_vars > 0 & !is.na(feature_vars)
  features_filtered <- features_data[, valid_features]
  
  cat(paste("   📊 Features:", ncol(features_data), "→", ncol(features_filtered), "(após filtragem)\n"))
  
  if(ncol(features_filtered) < 3) {
    cat(paste("   ⚠ Poucas features válidas para PCA\n\n"))
    return(NULL)
  }
  
  # 3. CLR transformation (adicionar pseudocount se necessário)
  # Adicionar pequeno valor para evitar zeros
  features_clr_input <- features_filtered + 1e-6
  
  # Aplicar CLR transformation
  features_clr <- compositions::clr(features_clr_input)
  
  # 4. Realizar PCA
  pca_result <- prcomp(features_clr, center = TRUE, scale. = TRUE)
  
  # Calcular variância explicada
  var_exp <- summary(pca_result)$importance["Proportion of Variance", 1:3] * 100
  
  # 5. Preparar dados para plot
  pca_df <- as.data.frame(pca_result$x[, 1:3]) %>%
    rownames_to_column("Row_Index") %>%
    mutate(Sample_ID = metadata_subset$Sample_ID,
           Group = metadata_subset[[group_var]],
           age = metadata_subset$age) %>%
    select(Sample_ID, Group, age, PC1, PC2, PC3)
  
  # 6. Criar gráfico PC1 vs PC2
  p_pca_12 <- ggplot(pca_df,
                     aes(x = PC1, y = PC2,
                         color = Group,
                         label = Sample_ID)) +
    geom_point(size = 3.5, alpha = 0.85) +
    geom_text_repel(
      size = 2.5,
      max.overlaps = 30,
      show.legend = FALSE,
      segment.color = "grey60"
    ) +
    stat_ellipse(
      aes(group = Group),
      type = "norm",
      linetype = "dashed",
      alpha = 0.5,
      linewidth = 0.7
    ) +
    scale_color_manual(values = c("control" = "#E31A1C", "glycodex" = "#1F78B4")) +
    labs(
      title = paste("PCA — Sample Swap Check:", dataset_name, "(PC1 vs PC2)"),
      subtitle = "Control vs Glycodex · CLR transformation",
      x = sprintf("PC1 (%.1f%%)", var_exp[1]),
      y = sprintf("PC2 (%.1f%%)", var_exp[2]),
      color = "Treatment"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey40"),
      legend.position = "right"
    )
  
  # 7. Criar gráfico PC1 vs PC3
  p_pca_13 <- ggplot(pca_df,
                     aes(x = PC1, y = PC3,
                         color = Group,
                         label = Sample_ID)) +
    geom_point(size = 3.5, alpha = 0.85) +
    geom_text_repel(
      size = 2.5,
      max.overlaps = 30,
      show.legend = FALSE,
      segment.color = "grey60"
    ) +
    stat_ellipse(
      aes(group = Group),
      type = "norm",
      linetype = "dashed",
      alpha = 0.5,
      linewidth = 0.7
    ) +
    scale_color_manual(values = c("control" = "#E31A1C", "glycodex" = "#1F78B4")) +
    labs(
      title = paste("PCA — Sample Swap Check:", dataset_name, "(PC1 vs PC3)"),
      subtitle = "Control vs Glycodex · CLR transformation",
      x = sprintf("PC1 (%.1f%%)", var_exp[1]),
      y = sprintf("PC3 (%.1f%%)", var_exp[3]),
      color = "Treatment"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey40"),
      legend.position = "right"
    )
  
  cat(paste("   ✅ PCA concluído - Variância explicada PC1-3:", 
            paste(sprintf("%.1f%%", var_exp), collapse = ", "), "\n"))
  
  # Distribuição por grupo
  group_counts <- table(pca_df$Group)
  cat(paste("   📊 Distribuição:", paste(names(group_counts), group_counts, sep = "=", collapse = ", "), "\n\n"))
  
  return(list(
    pca_result = pca_result,
    pca_data = pca_df,
    variance_explained = var_exp,
    plot_pc12 = p_pca_12,
    plot_pc13 = p_pca_13,
    dataset_name = dataset_name,
    group_distribution = group_counts
  ))
}

# Lista dos datasets AGRUPADOS para análise PCA
datasets_for_pca <- list(
  "KO" = df_ko_grouped,
  "AMR" = df_AMR_grouped,
  "VF" = df_VF_grouped,
  "CAZy" = df_CAZy_grouped,
  "EC" = df_EC_grouped,
  "COGs" = df_COGs_grouped
)

# Verificar disponibilidade dos datasets agrupados
cat("🎯 Verificando datasets agrupados disponíveis:\n")
for(name in names(datasets_for_pca)) {
  dataset <- datasets_for_pca[[name]]
  if(!is.null(dataset) && "treatment" %in% names(dataset)) {
    dims <- dim(dataset)
    groups <- table(dataset$treatment)
    cat(paste("✅", name, ":", dims[1], "amostras |", 
              paste(names(groups), groups, sep = "=", collapse = ", "), "\n"))
  } else {
    cat(paste("❌", name, ": não disponível\n"))
  }
}
cat("\n")

# Executar PCA para cada dataset agrupado
pca_results_grouped <- list()

for(dataset_name in names(datasets_for_pca)) {
  pca_results_grouped[[dataset_name]] <- perform_functional_pca(
    datasets_for_pca[[dataset_name]], 
    dataset_name
  )
}

# Exibir todos os gráficos PC1 vs PC2
cat("📊 EXIBINDO GRÁFICOS PC1 vs PC2 (CONTROL vs GLYCODEX):\n")
cat(rep("=", 50), "\n")

for(dataset_name in names(pca_results_grouped)) {
  if(!is.null(pca_results_grouped[[dataset_name]])) {
    print(pca_results_grouped[[dataset_name]]$plot_pc12)
    cat("\n")
  }
}


# Resumo da variância explicada para grupos selecionados
cat("📈 RESUMO DA VARIÂNCIA EXPLICADA (CONTROL vs GLYCODEX):\n")
cat(rep("=", 60), "\n")

variance_summary_grouped <- data.frame(
  Dataset = character(),
  PC1_percent = numeric(),
  PC2_percent = numeric(),
  PC3_percent = numeric(),
  PC123_total = numeric(),
  Control_n = integer(),
  Glycodex_n = integer(),
  stringsAsFactors = FALSE
)

for(dataset_name in names(pca_results_grouped)) {
  if(!is.null(pca_results_grouped[[dataset_name]])) {
    var_exp <- pca_results_grouped[[dataset_name]]$variance_explained
    group_dist <- pca_results_grouped[[dataset_name]]$group_distribution
    
    variance_summary_grouped <- rbind(variance_summary_grouped, data.frame(
      Dataset = dataset_name,
      PC1_percent = round(var_exp[1], 1),
      PC2_percent = round(var_exp[2], 1),
      PC3_percent = round(var_exp[3], 1),
      PC123_total = round(sum(var_exp[1:3]), 1),
      Control_n = group_dist["control"],
      Glycodex_n = group_dist["glycodex"]
    ))
  }
}

print(variance_summary_grouped)

# Análise de separação entre grupos
cat("\n🔍 ANÁLISE DE SEPARAÇÃO ENTRE GRUPOS:\n")
cat(rep("=", 50), "\n")

for(dataset_name in names(pca_results_grouped)) {
  if(!is.null(pca_results_grouped[[dataset_name]])) {
    pca_data <- pca_results_grouped[[dataset_name]]$pca_data
    
    # Calcular centróides dos grupos
    centroids <- pca_data %>%
      group_by(Group) %>%
      summarise(
        PC1_mean = mean(PC1),
        PC2_mean = mean(PC2),
        PC3_mean = mean(PC3),
        .groups = 'drop'
      )
    
    # Distância euclidiana entre centróides (PC1 vs PC2)
    if(nrow(centroids) == 2) {
      dist_pc12 <- sqrt((centroids$PC1_mean[1] - centroids$PC1_mean[2])^2 + 
                          (centroids$PC2_mean[1] - centroids$PC2_mean[2])^2)
      
      cat(paste(dataset_name, "- Distância entre centróides (PC1-PC2):", 
                round(dist_pc12, 2), "\n"))
    }
  }
}



for(dataset_name in names(pca_results_grouped)) {
  if(!is.null(pca_results_grouped[[dataset_name]])) {
    pca_data <- pca_results_grouped[[dataset_name]]$pca_data
    
    # Calcular distâncias do centróide por grupo
    outliers_detected <- pca_data %>%
      group_by(Group) %>%
      mutate(
        PC1_center = mean(PC1),
        PC2_center = mean(PC2),
        distance_from_center = sqrt((PC1 - PC1_center)^2 + (PC2 - PC2_center)^2),
        is_outlier = distance_from_center > (mean(distance_from_center) + 2*sd(distance_from_center))
      ) %>%
      filter(is_outlier) %>%
      select(Sample_ID, Group, distance_from_center)
    
    if(nrow(outliers_detected) > 0) {
      cat(paste(dataset_name, "- Possíveis outliers:\n"))
      print(outliers_detected)
      cat("\n")
    } else {
      cat(paste(dataset_name, "- Nenhum outlier detectado\n"))
    }
  }
}











































