###############################################################################
#                    ╔═══════════════════════════════════════╗
#                    ║   CONFIGURAÇÃO DO USUÁRIO - INÍCIO   ║
#                    ╚═══════════════════════════════════════╝
###############################################################################

# 1. DIRETÓRIO DOS DADOS
data_dir <- "C:/Users/dti-/Documents/ANH/Simphyome_Analysis/Hav1801"

# 2. COLUNA DE IDENTIFICAÇÃO DAS AMOSTRAS (no metadata)
sample_id_col_name <- "Sample_ID"

# 3. COLUNA DE GRUPOS E GRUPOS SELECIONADOS
group_col_name   <- "group"
selected_groups  <- c("control", "glycodex")

# 4. COLUNA DE IDADE E IDADES SELECIONADAS
#    - Se selected_ages = NULL -> mantém todas as idades
#    - Se selected_ages = c("P10", "P28") -> filtra apenas essas idades
age_col_name   <- "age"
selected_ages  <- NULL

# 5. CORES PARA OS GRUPOS NOS GRÁFICOS
group_colors <- c("#E31A1C", "#1F78B4")
names(group_colors) <- selected_groups

###############################################################################
#                    ╔═══════════════════════════════════════╗
#                    ║    CONFIGURAÇÃO DO USUÁRIO - FIM      ║
#                    ╚═══════════════════════════════════════╝
#  NAO ALTERE NADA ABAIXO DESTA LINHA
###############################################################################

# ═══════════════════════════════════════════════════════════════════════════════
# CARREGAMENTO DE BIBLIOTECAS
# ═══════════════════════════════════════════════════════════════════════════════

library(phyloseq)
library(DESeq2)
library(edgeR)
library(compositions)
library(zCompositions)
library(dplyr)
library(tibble)
library(readr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(corrplot)
library(VennDiagram)
library(readxl)
library(vegan)
library(ape)
library(ggpubr)
library(patchwork)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(openxlsx)
library(tidyr)

# Resolver conflitos de namespace
select    <- dplyr::select
filter    <- dplyr::filter
rename    <- dplyr::rename
mutate    <- dplyr::mutate
arrange   <- dplyr::arrange
summarise <- dplyr::summarise

# ═══════════════════════════════════════════════════════════════════════════════
# FASE 1 — IMPORTAÇÃO DOS DADOS
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n", rep("=", 60), "\n")
cat("=== FASE 1: IMPORTACAO DOS DADOS ===\n")
cat(rep("=", 60), "\n")
cat("Diretorio:", data_dir, "\n\n")

# Dicionário de padrões de arquivos
patterns <- c(
  "AMR"   = "AMR_counts",
  "CAZy"  = "CAZy_counts",
  "COGs"  = "COGs_count",
  "EC"    = "EC_count",
  "KO"    = "KO_counts",
  "VF"    = "VF_counts",
  "meta"  = "meta"
)

# Objeto de armazenamento
data_storage <- list()

# Busca e importação
for (key in names(patterns)) {
  file_found <- list.files(path = data_dir, pattern = patterns[key], full.names = TRUE)
  
  if (length(file_found) == 0) {
    warning(paste("No files found for:", key))
    next
  } else if (length(file_found) > 1) {
    message(paste("Multiple files found for:", key, "- using:", basename(file_found[1])))
    file_found <- file_found[1]
  }
  
  if (grepl("\\.xlsx$", file_found)) {
    data_storage[[key]] <- read_excel(file_found)
  } else {
    data_storage[[key]] <- read_delim(file_found, delim = "\t", show_col_types = FALSE)
  }
  
  message(paste("Imported:", key, "-", basename(file_found)))
}

# Atribuir aos dataframes individuais
df_ko   <- data_storage$KO
df_AMR  <- data_storage$AMR
df_VF   <- data_storage$VF
df_CAZy <- data_storage$CAZy
df_EC   <- data_storage$EC
df_COGs <- data_storage$COGs
df_meta <- data_storage$meta

# Padronizar nome da coluna de amostras no metadata
possible_id_cols <- c("sample", "Sample", "sample_id", "SampleID", "sample_ID", "Sample_ID")
id_col_found <- FALSE
for (col in possible_id_cols) {
  if (col %in% colnames(df_meta) && col != sample_id_col_name) {
    colnames(df_meta)[colnames(df_meta) == col] <- sample_id_col_name
    cat(paste("Coluna do metadata renomeada:", col, "->", sample_id_col_name, "\n"))
    id_col_found <- TRUE
    break
  } else if (col %in% colnames(df_meta) && col == sample_id_col_name) {
    cat(paste("Coluna", sample_id_col_name, "ja existe no metadata\n"))
    id_col_found <- TRUE
    break
  }
}

if (!id_col_found) {
  cat("AVISO: Nenhuma coluna de amostras reconhecida no metadata!\n")
  cat("Colunas disponiveis:", paste(colnames(df_meta), collapse = ", "), "\n")
}

# Verificação dos formatos originais
cat("\n=== VERIFICACAO DOS FORMATOS ORIGINAIS ===\n")
datasets <- list("KO" = df_ko, "AMR" = df_AMR, "VF" = df_VF,
                 "CAZy" = df_CAZy, "EC" = df_EC, "COGs" = df_COGs)

for (name in names(datasets)) {
  if (!is.null(datasets[[name]])) {
    dims <- dim(datasets[[name]])
    cat(paste(name, ":", dims[1], "linhas x", dims[2], "colunas\n"))
    cat("  Primeiras colunas:", paste(names(datasets[[name]])[1:min(3, ncol(datasets[[name]]))], collapse = ", "), "\n")
  }
}

# ═══════════════════════════════════════════════════════════════════════════════
# FASE 2 — TRANSPOSIÇÃO E PADRONIZAÇÃO
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n", rep("=", 60), "\n")
cat("=== FASE 2: TRANSPOSICAO E PADRONIZACAO ===\n")
cat(rep("=", 60), "\n")

# Função para transpor datasets funcionais
transpose_functional_data <- function(data, feature_id_name) {
  feature_ids   <- data[[1]]
  sample_names  <- names(data)[-1]
  count_matrix  <- data[, -1]
  
  transposed_matrix <- t(count_matrix)
  transposed_df     <- as.data.frame(transposed_matrix)
  names(transposed_df) <- feature_ids
  
  transposed_df[[sample_id_col_name]] <- sample_names
  
  transposed_df <- transposed_df[, c(sample_id_col_name,
                                     setdiff(names(transposed_df), sample_id_col_name))]
  
  return(transposed_df)
}

# Datasets que precisam de transposição
datasets_to_transpose <- list(
  "df_ko"   = list(data = df_ko,   label = "ko"),
  "df_CAZy" = list(data = df_CAZy, label = "cazy"),
  "df_COGs" = list(data = df_COGs, label = "cog"),
  "df_EC"   = list(data = df_EC,   label = "ec")
)

for (name in names(datasets_to_transpose)) {
  item <- datasets_to_transpose[[name]]
  if (!is.null(item$data)) {
    transposed <- transpose_functional_data(item$data, item$label)
    assign(name, transposed, envir = .GlobalEnv)
    cat(paste("OK", name, "transposto:", nrow(transposed), "x", ncol(transposed), "\n"))
  } else {
    cat(paste("AVISO:", name, "nao encontrado\n"))
  }
}

# Atualizar referências locais após transposição
df_ko   <- get("df_ko",   envir = .GlobalEnv)
df_CAZy <- get("df_CAZy", envir = .GlobalEnv)
df_COGs <- get("df_COGs", envir = .GlobalEnv)
df_EC   <- get("df_EC",   envir = .GlobalEnv)

# ═══════════════════════════════════════════════════════════════════════════════
# FASE 2B — PIVOTAR DATASETS EM FORMATO LONGO (AMR, VF)
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n", rep("=", 60), "\n")
cat("=== FASE 2B: PIVOTAR AMR E VF (FORMATO LONGO -> LARGO) ===\n")
cat(rep("=", 60), "\n\n")

# Função genérica para detectar e pivotar datasets longos
pivot_long_to_wide <- function(dataset, dataset_name, id_col) {
  
  if (is.null(dataset)) return(NULL)
  
  # Detectar se é formato longo: Sample_ID não é único
  if (id_col %in% names(dataset)) {
    n_rows   <- nrow(dataset)
    n_unique <- length(unique(dataset[[id_col]]))
    
    if (n_unique == n_rows) {
      cat(paste("  ", dataset_name, ": formato largo OK (", n_rows, "linhas,", n_unique, "IDs unicos)\n"))
      return(dataset)  # já está em formato largo
    }
    
    cat(paste("  ", dataset_name, ": FORMATO LONGO detectado (", n_rows, "linhas,", n_unique, "IDs unicos)\n"))
    cat(paste("    Colunas:", paste(names(dataset), collapse = ", "), "\n"))
    
  } else {
    # Se não tem a coluna de ID, verificar primeira coluna
    first_col <- names(dataset)[1]
    n_rows   <- nrow(dataset)
    n_unique <- length(unique(dataset[[first_col]]))
    
    if (n_unique == n_rows) {
      cat(paste("  ", dataset_name, ": formato largo OK\n"))
      return(dataset)
    }
    
    cat(paste("  ", dataset_name, ": FORMATO LONGO detectado pela coluna", first_col, "\n"))
  }
  
  # ── Identificar a coluna de feature (gene/template ID) ──────────────────
  # Procurar pela primeira coluna character que não é metadata
  possible_meta <- c(id_col, "group", "age", "treatment", "study", 
                     "tissue", "age_group", "sample")
  
  all_cols    <- names(dataset)
  non_meta    <- setdiff(all_cols, possible_meta)
  
  # Encontrar coluna de feature ID (primeira coluna character não-metadata)
  feature_col <- NULL
  for (col in non_meta) {
    if (is.character(dataset[[col]]) || is.factor(dataset[[col]])) {
      feature_col <- col
      break
    }
  }
  
  if (is.null(feature_col)) {
    cat(paste("    AVISO: Nao encontrou coluna de feature ID. Usando primeira coluna non-meta:", non_meta[1], "\n"))
    feature_col <- non_meta[1]
  }
  
  cat(paste("    Feature ID column:", feature_col, "\n"))
  cat(paste("    Features unicas:", length(unique(dataset[[feature_col]])), "\n"))
  
  # ── Identificar coluna de valor (score/count) ──────────────────────────
  # Procurar colunas numéricas que façam sentido como contagens
  numeric_cols <- non_meta[sapply(dataset[, non_meta, drop = FALSE], is.numeric)]
  
  # Prioridade: procurar por nomes comuns de contagem
  count_candidates <- c("Depth", "depth", "Hits", "hits", "Count", "count",
                        "Counts", "counts", "Abundance", "abundance",
                        "Score", "score", "Mapping_Length", "mapping_length")
  
  value_col <- NULL
  for (candidate in count_candidates) {
    if (candidate %in% numeric_cols) {
      value_col <- candidate
      break
    }
  }
  
  # Se não encontrou, contar ocorrências (presença)
  use_count <- FALSE
  if (is.null(value_col)) {
    cat("    Nenhuma coluna de valor padrao encontrada. Usando contagem de ocorrencias.\n")
    use_count <- TRUE
  } else {
    cat(paste("    Value column:", value_col, "\n"))
  }
  
  # ── Pivotar ─────────────────────────────────────────────────────────────
  
  if (use_count) {
    # Contar ocorrências de cada feature por amostra
    pivot_data <- dataset[, c(id_col, feature_col), drop = FALSE]
    pivot_data$count_val <- 1
    
    # Agregar: somar contagens
    agg_formula <- as.formula(paste("count_val ~", id_col, "+", feature_col))
    agg_data <- aggregate(agg_formula, data = pivot_data, FUN = sum)
    
    # Pivotar
    wide_df <- reshape(agg_data,
                       idvar     = id_col,
                       timevar   = feature_col,
                       direction = "wide",
                       v.names   = "count_val")
    
    # Limpar nomes das colunas (remover prefixo "count_val.")
    names(wide_df) <- gsub("^count_val\\.", "", names(wide_df))
    
  } else {
    # Usar a coluna de valor identificada
    pivot_data <- dataset[, c(id_col, feature_col, value_col), drop = FALSE]
    
    # Agregar: somar valores duplicados
    agg_formula <- as.formula(paste(value_col, "~", id_col, "+", feature_col))
    agg_data <- aggregate(agg_formula, data = pivot_data, FUN = sum)
    
    # Pivotar
    wide_df <- reshape(agg_data,
                       idvar     = id_col,
                       timevar   = feature_col,
                       direction = "wide",
                       v.names   = value_col)
    
    # Limpar nomes das colunas
    prefix <- paste0(value_col, ".")
    names(wide_df) <- gsub(paste0("^", gsub("([.])", "\\\\\\1", prefix)), "", names(wide_df))
  }
  
  # Substituir NAs por 0
  wide_df[is.na(wide_df)] <- 0
  
  # Remover atributo reshape
  attr(wide_df, "reshapeWide") <- NULL
  rownames(wide_df) <- NULL
  
  cat(paste("    Resultado:", nrow(wide_df), "amostras x", ncol(wide_df), "colunas\n"))
  cat(paste("    Primeiras colunas:", paste(names(wide_df)[1:min(5, ncol(wide_df))], collapse = ", "), "\n\n"))
  
  return(wide_df)
}

# ── Aplicar pivotação nos datasets que precisam ──────────────────────────────

# Verificar e pivotar df_AMR
df_AMR <- pivot_long_to_wide(df_AMR, "AMR", sample_id_col_name)

# Verificar e pivotar df_VF
df_VF <- pivot_long_to_wide(df_VF, "VF", sample_id_col_name)

# Verificar os outros (devem retornar sem mudanças)
df_ko   <- pivot_long_to_wide(df_ko,   "KO",   sample_id_col_name)
df_CAZy <- pivot_long_to_wide(df_CAZy, "CAZy", sample_id_col_name)
df_EC   <- pivot_long_to_wide(df_EC,   "EC",   sample_id_col_name)
df_COGs <- pivot_long_to_wide(df_COGs, "COGs", sample_id_col_name)

# ── Verificação final ─────────────────────────────────────────────────────────

cat("=== VERIFICACAO FINAL POS-PIVOT ===\n")
check_all <- list("KO" = df_ko, "AMR" = df_AMR, "VF" = df_VF,
                  "CAZy" = df_CAZy, "EC" = df_EC, "COGs" = df_COGs)

for (nm in names(check_all)) {
  df <- check_all[[nm]]
  if (!is.null(df)) {
    n_rows   <- nrow(df)
    n_unique <- length(unique(df[[sample_id_col_name]]))
    is_ok    <- n_rows == n_unique
    cat(paste(ifelse(is_ok, "OK", "ERRO"), nm, ":", 
              n_rows, "linhas |", n_unique, "IDs unicos |",
              ncol(df), "colunas\n"))
  }
}

cat("\nFASE 2B CONCLUIDA!\n")
# ═══════════════════════════════════════════════════════════════════════════════
# FASE 3 — PADRONIZAÇÃO DOS NOMES DAS AMOSTRAS
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n", rep("=", 60), "\n")
cat("=== FASE 3: PADRONIZACAO DOS NOMES DAS AMOSTRAS ===\n")
cat(rep("=", 60), "\n")

standardize_sample_names <- function(sample_names) {
  sample_names %>%
    str_replace_all("\\.", "-") %>%
    str_replace("_$", "")
}

# Padronizar nomes no metadata
if (!is.null(df_meta) && sample_id_col_name %in% names(df_meta)) {
  cat("Antes:", paste(head(df_meta[[sample_id_col_name]], 3), collapse = ", "), "\n")
  df_meta[[sample_id_col_name]] <- standardize_sample_names(df_meta[[sample_id_col_name]])
  cat("Depois:", paste(head(df_meta[[sample_id_col_name]], 3), collapse = ", "), "\n")
}

# Verificação de consistência dos Sample_IDs
cat("\nVerificacao de consistencia dos Sample_IDs:\n")

all_datasets <- list(
  "KO" = df_ko, "AMR" = df_AMR, "VF" = df_VF,
  "CAZy" = df_CAZy, "EC" = df_EC, "COGs" = df_COGs
)

all_sample_ids <- list()
for (name in names(all_datasets)) {
  if (!is.null(all_datasets[[name]]) && sample_id_col_name %in% names(all_datasets[[name]])) {
    all_sample_ids[[name]] <- all_datasets[[name]][[sample_id_col_name]]
  }
}
all_sample_ids[["META"]] <- df_meta[[sample_id_col_name]]

if (length(all_sample_ids) > 1) {
  first_ids <- sort(all_sample_ids[[1]])
  consistent <- TRUE
  for (i in 2:length(all_sample_ids)) {
    if (!identical(first_ids, sort(all_sample_ids[[i]]))) {
      consistent <- FALSE
      cat(paste("AVISO: Inconsistencia:", names(all_sample_ids)[1], "vs", names(all_sample_ids)[i], "\n"))
      diff_1 <- setdiff(all_sample_ids[[1]], all_sample_ids[[i]])
      diff_2 <- setdiff(all_sample_ids[[i]], all_sample_ids[[1]])
      if (length(diff_1) > 0) cat(paste("  Em", names(all_sample_ids)[1], "mas nao em", names(all_sample_ids)[i], ":", paste(head(diff_1, 5), collapse = ", "), "\n"))
      if (length(diff_2) > 0) cat(paste("  Em", names(all_sample_ids)[i], "mas nao em", names(all_sample_ids)[1], ":", paste(head(diff_2, 5), collapse = ", "), "\n"))
    }
  }
  if (consistent) cat("OK: Todos os Sample_IDs consistentes!\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# FASE 4 — FILTRAGEM DE AMOSTRAS PELO METADATA
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n", rep("=", 60), "\n")
cat("=== FASE 4: FILTRAGEM DE AMOSTRAS PELO METADATA ===\n")
cat(rep("=", 60), "\n")

reference_sample_ids <- df_meta[[sample_id_col_name]]
cat("Amostras de referencia no metadata:", length(reference_sample_ids), "\n\n")

filter_samples_by_metadata <- function(dataset, dataset_name, reference_ids) {
  if (is.null(dataset) || !sample_id_col_name %in% names(dataset)) {
    cat(paste("AVISO:", dataset_name, "- sem coluna", sample_id_col_name, "\n"))
    return(NULL)
  }
  
  original_n  <- nrow(dataset)
  filtered_df <- dataset[dataset[[sample_id_col_name]] %in% reference_ids, ]
  removed_n   <- original_n - nrow(filtered_df)
  
  cat(paste("  ", dataset_name, ":", original_n, "->", nrow(filtered_df),
            "(removidas:", removed_n, ")\n"))
  
  return(filtered_df)
}

df_ko   <- filter_samples_by_metadata(df_ko,   "df_ko",   reference_sample_ids)
df_AMR  <- filter_samples_by_metadata(df_AMR,  "df_AMR",  reference_sample_ids)
df_VF   <- filter_samples_by_metadata(df_VF,   "df_VF",   reference_sample_ids)
df_CAZy <- filter_samples_by_metadata(df_CAZy, "df_CAZy", reference_sample_ids)
df_EC   <- filter_samples_by_metadata(df_EC,   "df_EC",   reference_sample_ids)
df_COGs <- filter_samples_by_metadata(df_COGs, "df_COGs", reference_sample_ids)

# Verificação pós-filtragem
cat("\nVerificacao pos-filtragem:\n")
datasets_filtered <- list("META" = df_meta, "KO" = df_ko, "AMR" = df_AMR,
                          "VF" = df_VF, "CAZy" = df_CAZy, "EC" = df_EC, "COGs" = df_COGs)

sample_counts <- sapply(datasets_filtered, function(df) {
  if (!is.null(df)) nrow(df) else NA
})

for (nm in names(sample_counts)) {
  cat(paste("  ", nm, ":", sample_counts[nm], "amostras\n"))
}

if (length(unique(na.omit(sample_counts))) == 1) {
  cat("OK: Todos os datasets com", unique(na.omit(sample_counts)), "amostras\n")
} else {
  cat("AVISO: Numeros inconsistentes de amostras!\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# FASE 5 — JOIN DAS COLUNAS DO METADATA (group + age)
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n", rep("=", 60), "\n")
cat("=== FASE 5: JOIN DAS COLUNAS DO METADATA ===\n")
cat(rep("=", 60), "\n")

# Verificar que as colunas configuradas existem no metadata
required_meta_cols <- c(sample_id_col_name, group_col_name, age_col_name)
missing_cols <- setdiff(required_meta_cols, names(df_meta))

if (length(missing_cols) > 0) {
  stop(paste("ERRO: Colunas nao encontradas no metadata:", paste(missing_cols, collapse = ", "),
             "\n   Colunas disponiveis:", paste(names(df_meta), collapse = ", ")))
}

# Subset do metadata com as colunas de interesse
metadata_subset <- df_meta[, required_meta_cols, drop = FALSE]

cat("Colunas do metadata para join:", paste(required_meta_cols, collapse = ", "), "\n")
cat("Grupos disponiveis:", paste(unique(df_meta[[group_col_name]]), collapse = ", "), "\n")
cat("Idades disponiveis:", paste(sort(unique(df_meta[[age_col_name]])), collapse = ", "), "\n\n")

# Lista dos dataframes para join
dfs <- list(
  df_ko   = df_ko,
  df_AMR  = df_AMR,
  df_VF   = df_VF,
  df_CAZy = df_CAZy,
  df_EC   = df_EC,
  df_COGs = df_COGs
)

# Aplicar merge (base R para evitar conflitos)
dfs_joined <- lapply(names(dfs), function(name) {
  df <- dfs[[name]]
  if (is.null(df)) {
    cat(paste("AVISO:", name, "- nao encontrado\n"))
    return(NULL)
  }
  
  original_ncol <- ncol(df)
  
  df_with_meta <- merge(df, metadata_subset, by = sample_id_col_name, all.x = TRUE)
  
  # Reordenar: Sample_ID, group, age, ... features ...
  cols_order <- c(sample_id_col_name, group_col_name, age_col_name,
                  setdiff(names(df_with_meta), c(sample_id_col_name, group_col_name, age_col_name)))
  df_with_meta <- df_with_meta[, cols_order]
  
  na_group <- sum(is.na(df_with_meta[[group_col_name]]))
  na_age   <- sum(is.na(df_with_meta[[age_col_name]]))
  
  cat(paste("OK", name, ":", original_ncol, "->", ncol(df_with_meta), "colunas",
            "| NAs group:", na_group, "| NAs age:", na_age, "\n"))
  
  return(df_with_meta)
})

names(dfs_joined) <- names(dfs)
list2env(dfs_joined, envir = .GlobalEnv)

# Atualizar referências locais
df_ko   <- dfs_joined$df_ko
df_AMR  <- dfs_joined$df_AMR
df_VF   <- dfs_joined$df_VF
df_CAZy <- dfs_joined$df_CAZy
df_EC   <- dfs_joined$df_EC
df_COGs <- dfs_joined$df_COGs

# ═══════════════════════════════════════════════════════════════════════════════
# FASE 6 — FILTRO POR GRUPOS E IDADES SELECIONADOS
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n", rep("=", 60), "\n")
cat("=== FASE 6: FILTRO POR GRUPOS E IDADES ===\n")
cat(rep("=", 60), "\n")

cat("Grupos selecionados:", paste(selected_groups, collapse = ", "), "\n")
cat("Idades selecionadas:", ifelse(is.null(selected_ages), "TODAS", paste(selected_ages, collapse = ", ")), "\n\n")

# Distribuição atual
cat("Distribuicao dos grupos no metadata:\n")
print(table(df_meta[[group_col_name]]))
cat("\n")

if (!is.null(selected_ages)) {
  cat("Distribuicao das idades no metadata:\n")
  print(table(df_meta[[age_col_name]]))
  cat("\n")
}

# Lista dos dataframes para filtrar
dfs_to_filter <- list(
  df_ko   = df_ko,
  df_AMR  = df_AMR,
  df_VF   = df_VF,
  df_CAZy = df_CAZy,
  df_EC   = df_EC,
  df_COGs = df_COGs
)

# Aplicar filtros
dfs_grouped <- lapply(names(dfs_to_filter), function(name) {
  df <- dfs_to_filter[[name]]
  
  if (is.null(df) || !group_col_name %in% names(df)) {
    cat(paste("AVISO:", name, "- sem coluna", group_col_name, "\n"))
    return(NULL)
  }
  
  original_n <- nrow(df)
  
  # Filtro por grupos (base R)
  df_filtered <- df[df[[group_col_name]] %in% selected_groups, ]
  
  # Filtro por idades (se definido)
  if (!is.null(selected_ages)) {
    df_filtered <- df_filtered[df_filtered[[age_col_name]] %in% selected_ages, ]
  }
  
  new_n <- nrow(df_filtered)
  group_dist <- table(df_filtered[[group_col_name]])
  
  cat(paste("  ", name, ":", original_n, "->", new_n, "amostras |",
            paste(names(group_dist), group_dist, sep = "=", collapse = ", "), "\n"))
  
  if (!is.null(selected_ages)) {
    age_dist <- table(df_filtered[[age_col_name]])
    cat(paste("     Idades:", paste(names(age_dist), age_dist, sep = "=", collapse = ", "), "\n"))
  }
  
  return(df_filtered)
})

names(dfs_grouped) <- paste0(names(dfs_to_filter), "_grouped")
list2env(dfs_grouped, envir = .GlobalEnv)

# Atualizar referências locais
df_ko_grouped   <- dfs_grouped$df_ko_grouped
df_AMR_grouped  <- dfs_grouped$df_AMR_grouped
df_VF_grouped   <- dfs_grouped$df_VF_grouped
df_CAZy_grouped <- dfs_grouped$df_CAZy_grouped
df_EC_grouped   <- dfs_grouped$df_EC_grouped
df_COGs_grouped <- dfs_grouped$df_COGs_grouped

cat("\nFILTRO CONCLUIDO!\n")
cat("Datasets _grouped criados com:", paste(selected_groups, collapse = " & "), "\n")
if (!is.null(selected_ages)) {
  cat("Idades filtradas:", paste(selected_ages, collapse = ", "), "\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# FASE 7 — PCA PARA VERIFICAÇÃO DE TROCA DE AMOSTRAS
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n", rep("=", 60), "\n")
cat("=== FASE 7: PCA - SAMPLE SWAP CHECK ===\n")
cat(rep("=", 60), "\n")

# Colunas de metadata dentro dos datasets
meta_cols_internal <- c(sample_id_col_name, group_col_name, age_col_name)

# Função PCA parametrizada
perform_functional_pca <- function(dataset, dataset_name) {
  
  if (is.null(dataset) || !group_col_name %in% names(dataset)) {
    cat(paste("AVISO:", dataset_name, "- nao disponivel\n"))
    return(NULL)
  }
  
  cat(paste("PCA:", dataset_name, "...\n"))
  
  # Separar metadata das features
  metadata_cols_present <- intersect(meta_cols_internal, names(dataset))
  feature_cols  <- setdiff(names(dataset), meta_cols_internal)
  
  metadata_sub  <- dataset[, metadata_cols_present, drop = FALSE]
  features_data <- dataset[, feature_cols, drop = FALSE]
  
  # Converter para numérico
  features_data <- as.data.frame(lapply(features_data, function(x) {
    as.numeric(as.character(x))
  }))
  
  # Substituir NAs por 0
  features_data[is.na(features_data)] <- 0
  
  # Remover features com variância zero
  feature_vars     <- apply(features_data, 2, var, na.rm = TRUE)
  valid_features   <- !is.na(feature_vars) & feature_vars > 0
  features_filtered <- features_data[, valid_features, drop = FALSE]
  
  cat(paste("   Features:", ncol(features_data), "->", ncol(features_filtered), "\n"))
  
  if (ncol(features_filtered) < 3) {
    cat("   AVISO: Poucas features para PCA\n")
    return(NULL)
  }
  
  # CLR transformation
  features_clr <- compositions::clr(features_filtered + 1e-6)
  
  # PCA
  pca_result <- prcomp(features_clr, center = TRUE, scale. = TRUE)
  
  n_pcs <- min(3, ncol(pca_result$x))
  var_exp <- summary(pca_result)$importance["Proportion of Variance", 1:n_pcs] * 100
  
  # Dados para plot
  pca_df <- as.data.frame(pca_result$x[, 1:n_pcs, drop = FALSE])
  pca_df$Sample_ID <- metadata_sub[[sample_id_col_name]]
  pca_df$Group     <- metadata_sub[[group_col_name]]
  
  if (age_col_name %in% names(metadata_sub)) {
    pca_df$Age <- metadata_sub[[age_col_name]]
  }
  
  # Plot PC1 vs PC2
  p_pca_12 <- ggplot(pca_df,
                     aes(x = PC1, y = PC2, color = Group, label = Sample_ID)) +
    geom_point(size = 3.5, alpha = 0.85) +
    geom_text_repel(size = 2.5, max.overlaps = 30, show.legend = FALSE,
                    segment.color = "grey60") +
    stat_ellipse(aes(group = Group), type = "norm", linetype = "dashed",
                 alpha = 0.5, linewidth = 0.7) +
    scale_color_manual(values = group_colors) +
    labs(
      title    = paste("PCA - Sample Swap Check:", dataset_name, "(PC1 vs PC2)"),
      subtitle = paste(paste(selected_groups, collapse = " vs "), "- CLR transformation"),
      x = sprintf("PC1 (%.1f%%)", var_exp[1]),
      y = sprintf("PC2 (%.1f%%)", var_exp[2]),
      color = group_col_name
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(color = "grey40"),
          legend.position = "right")
  
  # Plot PC1 vs PC3 (se tiver 3+ PCs)
  p_pca_13 <- NULL
  if (n_pcs >= 3) {
    p_pca_13 <- ggplot(pca_df,
                       aes(x = PC1, y = PC3, color = Group, label = Sample_ID)) +
      geom_point(size = 3.5, alpha = 0.85) +
      geom_text_repel(size = 2.5, max.overlaps = 30, show.legend = FALSE,
                      segment.color = "grey60") +
      stat_ellipse(aes(group = Group), type = "norm", linetype = "dashed",
                   alpha = 0.5, linewidth = 0.7) +
      scale_color_manual(values = group_colors) +
      labs(
        title    = paste("PCA - Sample Swap Check:", dataset_name, "(PC1 vs PC3)"),
        subtitle = paste(paste(selected_groups, collapse = " vs "), "- CLR transformation"),
        x = sprintf("PC1 (%.1f%%)", var_exp[1]),
        y = sprintf("PC3 (%.1f%%)", var_exp[3]),
        color = group_col_name
      ) +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold"),
            plot.subtitle = element_text(color = "grey40"),
            legend.position = "right")
  }
  
  cat(paste("   OK Variancia PC1-", n_pcs, ":",
            paste(sprintf("%.1f%%", var_exp), collapse = ", "), "\n"))
  
  group_counts <- table(pca_df$Group)
  cat(paste("   ", paste(names(group_counts), group_counts, sep = "=", collapse = ", "), "\n\n"))
  
  return(list(
    pca_result         = pca_result,
    pca_data           = pca_df,
    variance_explained = var_exp,
    plot_pc12          = p_pca_12,
    plot_pc13          = p_pca_13,
    dataset_name       = dataset_name,
    group_distribution = group_counts
  ))
}

# Datasets agrupados para PCA
datasets_for_pca <- list(
  "KO"   = df_ko_grouped,
  "AMR"  = df_AMR_grouped,
  "VF"   = df_VF_grouped,
  "CAZy" = df_CAZy_grouped,
  "EC"   = df_EC_grouped,
  "COGs" = df_COGs_grouped
)

# Verificar disponibilidade
cat("Datasets disponiveis para PCA:\n")
for (name in names(datasets_for_pca)) {
  dataset <- datasets_for_pca[[name]]
  if (!is.null(dataset) && group_col_name %in% names(dataset)) {
    dims <- dim(dataset)
    groups <- table(dataset[[group_col_name]])
    cat(paste("  OK", name, ":", dims[1], "amostras |",
              paste(names(groups), groups, sep = "=", collapse = ", "), "\n"))
  } else {
    cat(paste("  AVISO:", name, "nao disponivel\n"))
  }
}
cat("\n")

# Executar PCA
pca_results_grouped <- list()
for (dataset_name in names(datasets_for_pca)) {
  pca_results_grouped[[dataset_name]] <- perform_functional_pca(
    datasets_for_pca[[dataset_name]], dataset_name
  )
}

# Exibir gráficos PC1 vs PC2
cat("=== GRAFICOS PC1 vs PC2 ===\n")
for (dataset_name in names(pca_results_grouped)) {
  if (!is.null(pca_results_grouped[[dataset_name]])) {
    print(pca_results_grouped[[dataset_name]]$plot_pc12)
  }
}

# Resumo de variância explicada
cat("\n=== RESUMO DA VARIANCIA EXPLICADA ===\n")
cat(rep("=", 60), "\n")

variance_summary <- do.call(rbind, lapply(names(pca_results_grouped), function(name) {
  res <- pca_results_grouped[[name]]
  if (!is.null(res)) {
    ve <- res$variance_explained
    n_pcs <- length(ve)
    data.frame(
      Dataset     = name,
      PC1         = round(ve[1], 1),
      PC2         = if (n_pcs >= 2) round(ve[2], 1) else NA,
      PC3         = if (n_pcs >= 3) round(ve[3], 1) else NA,
      PC123_total = round(sum(ve[1:min(3, n_pcs)]), 1),
      stringsAsFactors = FALSE
    )
  }
}))

# Adicionar contagens de grupo
for (g in selected_groups) {
  variance_summary[[paste0(g, "_n")]] <- sapply(variance_summary$Dataset, function(d) {
    gd <- pca_results_grouped[[d]]$group_distribution
    if (g %in% names(gd)) as.integer(gd[g]) else 0L
  })
}

print(variance_summary)

# Análise de separação entre grupos
cat("\n=== DISTANCIA ENTRE CENTROIDES (PC1-PC2) ===\n")

for (dataset_name in names(pca_results_grouped)) {
  res <- pca_results_grouped[[dataset_name]]
  if (!is.null(res)) {
    centroids <- aggregate(
      cbind(PC1, PC2) ~ Group,
      data = res$pca_data,
      FUN = mean
    )
    
    if (nrow(centroids) >= 2) {
      dist_pc12 <- sqrt((centroids$PC1[1] - centroids$PC1[2])^2 +
                          (centroids$PC2[1] - centroids$PC2[2])^2)
      cat(paste("  ", dataset_name, ":", round(dist_pc12, 2), "\n"))
    }
  }
}

# Detecção de outliers
cat("\n=== DETECCAO DE OUTLIERS (>2 SD do centroide) ===\n")

for (dataset_name in names(pca_results_grouped)) {
  res <- pca_results_grouped[[dataset_name]]
  if (!is.null(res)) {
    pca_data <- res$pca_data
    
    centroids <- aggregate(
      cbind(PC1, PC2) ~ Group,
      data = pca_data,
      FUN = mean
    )
    names(centroids)[2:3] <- c("PC1_center", "PC2_center")
    
    pca_with_center <- merge(pca_data, centroids, by = "Group")
    
    pca_with_center$distance <- sqrt(
      (pca_with_center$PC1 - pca_with_center$PC1_center)^2 +
        (pca_with_center$PC2 - pca_with_center$PC2_center)^2
    )
    
    outliers_list <- list()
    for (grp in unique(pca_with_center$Group)) {
      grp_data <- pca_with_center[pca_with_center$Group == grp, ]
      threshold <- mean(grp_data$distance) + 2 * sd(grp_data$distance)
      grp_outliers <- grp_data[grp_data$distance > threshold, ]
      if (nrow(grp_outliers) > 0) {
        outliers_list[[grp]] <- grp_outliers
      }
    }
    
    if (length(outliers_list) > 0) {
      all_outliers <- do.call(rbind, outliers_list)
      cat(paste("  ", dataset_name, "- Possiveis outliers:\n"))
      print(all_outliers[, c("Sample_ID", "Group", "distance")])
      cat("\n")
    } else {
      cat(paste("  ", dataset_name, "- Nenhum outlier\n"))
    }
  }
}

cat("\nFASE 7 CONCLUIDA!\n")

# ═══════════════════════════════════════════════════════════════════════════════
# ███████████████████████████████████████████████████████████████████████████████
# FASE 8 — ALPHA DIVERSITY
# ███████████████████████████████████████████████████████████████████████████████
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n", rep("=", 70), "\n")
cat("=== FASE 8: ALPHA DIVERSITY — TODOS OS DATASETS FUNCIONAIS ===\n")
cat(rep("=", 70), "\n\n")

# ── Setup ─────────────────────────────────────────────────────────────────────

# Lista central de datasets agrupados
datasets_grouped_list <- list(
  "KO"   = if (exists("df_ko_grouped"))   df_ko_grouped   else NULL,
  "AMR"  = if (exists("df_AMR_grouped"))  df_AMR_grouped  else NULL,
  "VF"   = if (exists("df_VF_grouped"))   df_VF_grouped   else NULL,
  "CAZy" = if (exists("df_CAZy_grouped")) df_CAZy_grouped else NULL,
  "EC"   = if (exists("df_EC_grouped"))   df_EC_grouped   else NULL,
  "COGs" = if (exists("df_COGs_grouped")) df_COGs_grouped else NULL
)
datasets_grouped_list <- datasets_grouped_list[!sapply(datasets_grouped_list, is.null)]

cat("Datasets disponiveis para analise:\n")
for (nm in names(datasets_grouped_list)) {
  dims <- dim(datasets_grouped_list[[nm]])
  cat(paste("  ", nm, ":", dims[1], "amostras x", dims[2], "colunas\n"))
}

# Diretório de output
output_root <- file.path(data_dir, "functional_diversity_results")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

# Paleta de cores expandível
n_groups <- length(selected_groups)
if (n_groups <= 2) {
  group_palette <- setNames(c("#E31A1C", "#1F78B4")[1:n_groups], selected_groups)
} else {
  group_palette <- setNames(
    RColorBrewer::brewer.pal(min(n_groups, 9), "Set1")[1:n_groups],
    selected_groups
  )
}

# Pares de comparação para stat_compare_means
if (n_groups == 2) {
  comparisons_list <- list(selected_groups)
} else {
  comparisons_list <- combn(selected_groups, 2, simplify = FALSE)
}

# ── Funções auxiliares ────────────────────────────────────────────────────────

extract_count_matrix <- function(dataset) {
  feature_cols <- setdiff(names(dataset), meta_cols_internal)
  mat <- as.data.frame(dataset[, feature_cols])
  mat <- as.data.frame(lapply(mat, function(x) as.numeric(as.character(x))))
  mat[is.na(mat)] <- 0
  rownames(mat) <- dataset[[sample_id_col_name]]
  as.matrix(mat)
}

extract_group_vector <- function(dataset) {
  grp <- dataset[[group_col_name]]
  names(grp) <- dataset[[sample_id_col_name]]
  grp
}

sig_stars <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "ns")))
}

save_gg <- function(plot_obj, filepath, width = 10, height = 7, dpi = 300) {
  ggsave(paste0(filepath, ".png"),
         plot   = plot_obj,
         width  = width,
         height = height,
         dpi    = dpi,
         bg     = "white")
  cat(paste("    Salvo:", basename(filepath), ".png\n"))
}

# ── Armazenamento de resultados ──────────────────────────────────────────────

alpha_results_all    <- list()
alpha_stats_all      <- list()
alpha_plots_all      <- list()
rarefaction_info_all <- list()

# ── Loop principal Alpha ─────────────────────────────────────────────────────

for (ds_name in names(datasets_grouped_list)) {
  
  cat(rep("-", 60), "\n")
  cat(paste("ALPHA DIVERSITY:", ds_name, "\n"))
  cat(rep("-", 60), "\n")
  
  dataset <- datasets_grouped_list[[ds_name]]
  
  # 8.1 Extrair matriz de contagens
  count_mat <- extract_count_matrix(dataset)
  group_vec <- extract_group_vector(dataset)
  
  # Arredondar para inteiros (necessário para rarefação)
  count_mat <- round(count_mat)
  count_mat[count_mat < 0] <- 0
  
  # Remover features com soma zero
  feat_sums <- colSums(count_mat)
  count_mat <- count_mat[, feat_sums > 0, drop = FALSE]
  
  cat(paste("  Matriz:", nrow(count_mat), "amostras x", ncol(count_mat), "features\n"))
  
  # 8.2 Profundidade de sequenciamento
  sample_depths <- rowSums(count_mat)
  rare_depth    <- min(sample_depths)
  
  cat(paste("  Profundidade min:", rare_depth,
            "| max:", max(sample_depths),
            "| mediana:", median(sample_depths), "\n"))
  
  rarefaction_info_all[[ds_name]] <- data.frame(
    Dataset      = ds_name,
    Min_Depth    = rare_depth,
    Max_Depth    = max(sample_depths),
    Median_Depth = median(sample_depths),
    N_Samples    = nrow(count_mat),
    N_Features   = ncol(count_mat),
    stringsAsFactors = FALSE
  )
  
  # Plot de profundidade
  depth_df <- data.frame(
    Sample = names(sample_depths),
    Depth  = as.numeric(sample_depths),
    Group  = group_vec[names(sample_depths)],
    stringsAsFactors = FALSE
  )
  
  p_depth <- ggplot(depth_df,
                    aes(x = reorder(Sample, Depth), y = Depth, fill = Group)) +
    geom_col(alpha = 0.85) +
    geom_hline(yintercept = rare_depth, color = "darkred",
               linetype = "dashed", linewidth = 0.8) +
    annotate("text", x = 1, y = rare_depth,
             label = sprintf("Rarefaction: %d", rare_depth),
             vjust = -0.5, hjust = 0, color = "darkred", size = 3.5) +
    scale_fill_manual(values = group_palette) +
    labs(
      title    = paste(ds_name, "- Sequencing Depth per Sample"),
      subtitle = "Red dashed line = rarefaction threshold (minimum depth)",
      x = "Sample", y = "Total Counts", fill = group_col_name
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title  = element_text(face = "bold"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
    )
  
  # 8.3 Rarefação
  if (rare_depth < 10) {
    cat("  AVISO: Profundidade muito baixa para rarefacao. Usando dados nao rarefeitos.\n")
    count_rare <- count_mat
  } else {
    set.seed(42)
    count_rare <- vegan::rrarefy(count_mat, sample = rare_depth)
    count_rare <- count_rare[, colSums(count_rare) > 0, drop = FALSE]
    cat(paste("  Apos rarefacao:", nrow(count_rare), "amostras x",
              ncol(count_rare), "features\n"))
  }
  
  # 8.4 Calcular métricas de alfa diversidade
  obs_richness <- rowSums(count_rare > 0)
  
  chao1_vals <- tryCatch({
    est <- t(vegan::estimateR(count_rare))
    est[, "S.chao1"]
  }, error = function(e) {
    cat("  AVISO: Chao1 falhou, usando Observed.\n")
    obs_richness
  })
  
  shannon_vals <- vegan::diversity(count_rare, index = "shannon")
  simpson_vals <- vegan::diversity(count_rare, index = "simpson")
  
  alpha_df <- data.frame(
    Sample   = rownames(count_rare),
    Group    = group_vec[rownames(count_rare)],
    Observed = as.numeric(obs_richness),
    Chao1    = as.numeric(chao1_vals),
    Shannon  = as.numeric(shannon_vals),
    Simpson  = as.numeric(simpson_vals),
    stringsAsFactors = FALSE
  )
  
  alpha_results_all[[ds_name]] <- alpha_df
  
  # 8.5 Formato longo para plots
  alpha_long <- tidyr::pivot_longer(
    alpha_df,
    cols      = c("Observed", "Chao1", "Shannon", "Simpson"),
    names_to  = "Metric",
    values_to = "Value"
  )
  alpha_long$Metric <- factor(alpha_long$Metric,
                              levels = c("Observed", "Chao1", "Shannon", "Simpson"))
  
  # 8.6 Plot facetado — 4 métricas
  p_alpha_facet <- ggplot(alpha_long,
                          aes(x = Group, y = Value,
                              fill = Group, color = Group)) +
    geom_violin(alpha = 0.35, linewidth = 0.6, trim = FALSE) +
    geom_boxplot(width = 0.18, alpha = 0.80,
                 outlier.shape = NA, color = "black") +
    geom_jitter(width = 0.08, size = 1.8, alpha = 0.70, show.legend = FALSE) +
    stat_compare_means(
      method       = "wilcox.test",
      comparisons  = comparisons_list,
      label        = "p.format",
      tip.length   = 0.01,
      bracket.size = 0.4,
      size         = 3.5
    ) +
    scale_fill_manual(values  = group_palette) +
    scale_color_manual(values = group_palette) +
    facet_wrap(~ Metric, scales = "free_y", nrow = 2) +
    labs(
      title    = paste("Alpha Diversity -", ds_name),
      subtitle = sprintf("Rarefied to %d counts/sample | Wilcoxon p-value",
                         rare_depth),
      x = NULL, y = "Diversity Value",
      fill = group_col_name, color = group_col_name
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title       = element_text(face = "bold"),
      plot.subtitle    = element_text(color = "grey40", size = 9.5),
      strip.text       = element_text(face = "bold", size = 12),
      strip.background = element_rect(fill = "grey92", color = NA),
      legend.position  = "bottom",
      axis.text.x      = element_text(angle = 20, hjust = 1)
    )
  
  # 8.7 Plots individuais 2x2
  plot_alpha_metric <- function(metric_name, y_label) {
    dat <- alpha_long[alpha_long$Metric == metric_name, ]
    
    ggplot(dat, aes(x = Group, y = Value, fill = Group, color = Group)) +
      geom_violin(alpha = 0.30, linewidth = 0.6, trim = FALSE) +
      geom_boxplot(width = 0.20, alpha = 0.80,
                   outlier.shape = NA, color = "black") +
      geom_jitter(width = 0.08, size = 2.2, alpha = 0.75, show.legend = FALSE) +
      stat_compare_means(
        method      = "wilcox.test",
        comparisons = comparisons_list,
        label       = "p.format",
        tip.length  = 0.01,
        size        = 3.8
      ) +
      scale_fill_manual(values  = group_palette) +
      scale_color_manual(values = group_palette) +
      labs(title = metric_name, x = NULL, y = y_label) +
      theme_bw(base_size = 13) +
      theme(
        plot.title      = element_text(face = "bold", hjust = 0.5),
        legend.position = "none",
        axis.text.x     = element_text(angle = 20, hjust = 1)
      )
  }
  
  p_obs <- plot_alpha_metric("Observed", "Observed Richness")
  p_ch  <- plot_alpha_metric("Chao1",    "Chao1 Estimated Richness")
  p_sha <- plot_alpha_metric("Shannon",  "Shannon Index")
  p_sim <- plot_alpha_metric("Simpson",  "Simpson Index")
  
  p_alpha_panel <- (p_obs | p_ch) / (p_sha | p_sim) +
    plot_annotation(
      title    = paste("Alpha Diversity -", ds_name, "- Individual Metrics"),
      subtitle = sprintf("Rarefied to %d counts/sample | Wilcoxon p-value",
                         rare_depth),
      theme = theme(
        plot.title    = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(color = "grey40", size = 10)
      )
    )
  
  # 8.8 Testes estatísticos
  metrics_list <- c("Observed", "Chao1", "Shannon", "Simpson")
  
  # Kruskal-Wallis
  kw_results <- do.call(rbind, lapply(metrics_list, function(m) {
    test <- kruskal.test(alpha_df[[m]] ~ factor(alpha_df$Group))
    data.frame(
      Dataset   = ds_name,
      Metric    = m,
      Test      = "Kruskal-Wallis",
      Statistic = round(test$statistic, 4),
      df        = test$parameter,
      p_value   = round(test$p.value, 4),
      Sig       = sig_stars(test$p.value),
      stringsAsFactors = FALSE
    )
  }))
  
  # Wilcoxon pairwise
  wx_results <- do.call(rbind, lapply(metrics_list, function(m) {
    do.call(rbind, lapply(comparisons_list, function(pair) {
      vals_a <- alpha_df[alpha_df$Group == pair[1], m]
      vals_b <- alpha_df[alpha_df$Group == pair[2], m]
      test   <- wilcox.test(vals_a, vals_b, exact = FALSE)
      data.frame(
        Dataset    = ds_name,
        Metric     = m,
        Test       = "Wilcoxon",
        Comparison = paste(pair[1], "vs", pair[2]),
        W          = round(test$statistic, 4),
        p_value    = round(test$p.value, 4),
        Sig        = sig_stars(test$p.value),
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  # Estatísticas descritivas
  desc_stats <- do.call(rbind, lapply(metrics_list, function(m) {
    do.call(rbind, lapply(selected_groups, function(g) {
      vals <- alpha_df[alpha_df$Group == g, m]
      data.frame(
        Dataset = ds_name,
        Group   = g,
        Metric  = m,
        N       = length(vals),
        Mean    = round(mean(vals, na.rm = TRUE), 4),
        Median  = round(median(vals, na.rm = TRUE), 4),
        SD      = round(sd(vals, na.rm = TRUE), 4),
        Min     = round(min(vals, na.rm = TRUE), 4),
        Max     = round(max(vals, na.rm = TRUE), 4),
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  alpha_stats_all[[ds_name]] <- list(
    kruskal_wallis = kw_results,
    wilcoxon       = wx_results,
    descriptive    = desc_stats
  )
  
  # 8.9 Armazenar plots
  alpha_plots_all[[ds_name]] <- list(
    depth    = p_depth,
    facet    = p_alpha_facet,
    panel    = p_alpha_panel,
    observed = p_obs,
    chao1    = p_ch,
    shannon  = p_sha,
    simpson  = p_sim
  )
  
  # 8.10 Imprimir resultados
  cat("\n  Kruskal-Wallis:\n")
  print(kw_results[, c("Metric", "Statistic", "p_value", "Sig")])
  cat("\n  Wilcoxon:\n")
  print(wx_results[, c("Metric", "Comparison", "W", "p_value", "Sig")])
  
  # Exibir plots
  print(p_depth)
  print(p_alpha_facet)
  print(p_alpha_panel)
  
  cat(paste("\n  OK:", ds_name, "alpha diversity completa!\n\n"))
}

# 8.11 Tabela resumo consolidada
cat("\n", rep("=", 70), "\n")
cat("=== RESUMO CONSOLIDADO - ALPHA DIVERSITY ===\n")
cat(rep("=", 70), "\n\n")

all_kw   <- do.call(rbind, lapply(alpha_stats_all, function(x) x$kruskal_wallis))
all_wx   <- do.call(rbind, lapply(alpha_stats_all, function(x) x$wilcoxon))
all_desc <- do.call(rbind, lapply(alpha_stats_all, function(x) x$descriptive))

cat("Kruskal-Wallis (todos os datasets):\n")
print(all_kw)
cat("\nWilcoxon (todos os datasets):\n")
print(all_wx)

# ═══════════════════════════════════════════════════════════════════════════════
# ███████████████████████████████████████████████████████████████████████████████
# FASE 9 — BETA DIVERSITY
# ███████████████████████████████████████████████████████████████████████████████
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n", rep("=", 70), "\n")
cat("=== FASE 9: BETA DIVERSITY — TODOS OS DATASETS FUNCIONAIS ===\n")
cat(rep("=", 70), "\n\n")

# Armazenamento
beta_stats_all  <- list()
beta_plots_all  <- list()

for (ds_name in names(datasets_grouped_list)) {
  
  cat(rep("-", 60), "\n")
  cat(paste("BETA DIVERSITY:", ds_name, "\n"))
  cat(rep("-", 60), "\n")
  
  dataset <- datasets_grouped_list[[ds_name]]
  
  # 9.1 Preparar dados
  count_mat <- extract_count_matrix(dataset)
  group_vec <- extract_group_vector(dataset)
  
  count_mat <- round(count_mat)
  count_mat[count_mat < 0] <- 0
  count_mat <- count_mat[, colSums(count_mat) > 0, drop = FALSE]
  
  # Rarefação para beta diversidade
  rare_depth <- min(rowSums(count_mat))
  
  if (rare_depth >= 10) {
    set.seed(42)
    count_rare <- vegan::rrarefy(count_mat, sample = rare_depth)
    count_rare <- count_rare[, colSums(count_rare) > 0, drop = FALSE]
  } else {
    count_rare <- count_mat
    cat("  AVISO: Profundidade baixa, usando dados nao rarefeitos.\n")
  }
  
  # Abundância relativa
  rel_mat <- sweep(count_rare, 1, rowSums(count_rare), "/")
  rel_mat[is.nan(rel_mat)] <- 0
  
  # Metadata alinhado
  meta_beta <- data.frame(
    Sample = rownames(count_rare),
    Group  = group_vec[rownames(count_rare)],
    stringsAsFactors = FALSE
  )
  
  cat(paste("  Matriz:", nrow(count_rare), "amostras x",
            ncol(count_rare), "features | Rarefied:", rare_depth, "\n"))
  
  # 9.2 Calcular matrizes de distância
  dist_bray <- vegdist(rel_mat, method = "bray")
  dist_jacc <- vegdist(count_rare, method = "jaccard", binary = TRUE)
  
  dist_list <- list(
    "Bray-Curtis" = dist_bray,
    "Jaccard"     = dist_jacc
  )
  
  # 9.3 PERMANOVA + ANOSIM
  cat("\n  Testes estatisticos:\n")
  
  stats_by_dist <- list()
  
  for (dist_name in names(dist_list)) {
    
    # PERMANOVA
    set.seed(42)
    perm <- adonis2(
      dist_list[[dist_name]] ~ Group,
      data         = meta_beta,
      permutations = 999,
      by           = "margin"
    )
    perm_R2 <- round(perm$R2[1], 3)
    perm_F  <- round(perm$F[1], 3)
    perm_p  <- perm$`Pr(>F)`[1]
    
    # ANOSIM
    set.seed(42)
    ano   <- anosim(dist_list[[dist_name]],
                    grouping     = meta_beta$Group,
                    permutations = 999)
    ano_R <- round(ano$statistic, 3)
    ano_p <- ano$signif
    
    annot_text <- sprintf(
      "PERMANOVA: R2=%.3f, F=%.3f, p=%.4f %s\nANOSIM: R=%.3f, p=%.4f %s",
      perm_R2, perm_F, perm_p, sig_stars(perm_p),
      ano_R, ano_p, sig_stars(ano_p)
    )
    
    cat(sprintf("    [%s] PERMANOVA R2=%.3f F=%.3f p=%.4f%s | ANOSIM R=%.3f p=%.4f%s\n",
                dist_name,
                perm_R2, perm_F, perm_p, sig_stars(perm_p),
                ano_R, ano_p, sig_stars(ano_p)))
    
    stats_by_dist[[dist_name]] <- list(
      permanova_R2  = perm_R2,
      permanova_F   = perm_F,
      permanova_p   = perm_p,
      permanova_sig = sig_stars(perm_p),
      anosim_R      = ano_R,
      anosim_p      = ano_p,
      anosim_sig    = sig_stars(ano_p),
      annot_text    = annot_text
    )
  }
  
  # 9.4 Betadisper (homogeneidade de dispersão)
  betadisp_results <- do.call(rbind, lapply(names(dist_list), function(dist_name) {
    beta_d <- betadisper(dist_list[[dist_name]], group = meta_beta$Group)
    set.seed(42)
    perm_d <- permutest(beta_d, permutations = 999)
    data.frame(
      Dataset  = ds_name,
      Distance = dist_name,
      F_stat   = round(perm_d$tab[1, "F"], 4),
      p_value  = round(perm_d$tab[1, "Pr(>F)"], 4),
      Result   = ifelse(perm_d$tab[1, "Pr(>F)"] < 0.05,
                        "Dispersao desigual - PERMANOVA pode refletir spread",
                        "Dispersao homogenea - PERMANOVA confiavel"),
      stringsAsFactors = FALSE
    )
  }))
  
  cat("\n  Betadisper:\n")
  print(betadisp_results[, c("Distance", "F_stat", "p_value", "Result")])
  
  # 9.5 PCoA
  pcoa_plots <- list()
  
  for (dist_name in names(dist_list)) {
    
    pcoa_res <- pcoa(dist_list[[dist_name]])
    eig      <- pcoa_res$values$Eigenvalues
    var_pcoa <- round(eig[eig > 0] / sum(eig[eig > 0]) * 100, 1)
    
    pcoa_df <- as.data.frame(pcoa_res$vectors[, 1:min(2, ncol(pcoa_res$vectors))])
    colnames(pcoa_df) <- c("Axis1", "Axis2")[1:ncol(pcoa_df)]
    pcoa_df$Sample <- rownames(pcoa_df)
    pcoa_df$Group  <- meta_beta$Group[match(pcoa_df$Sample, meta_beta$Sample)]
    
    annot <- stats_by_dist[[dist_name]]$annot_text
    
    p_pcoa <- ggplot(pcoa_df,
                     aes(x = Axis1, y = Axis2,
                         color = Group, label = Sample)) +
      geom_point(size = 3.5, alpha = 0.85) +
      geom_text_repel(size = 2.5, max.overlaps = 25,
                      show.legend = FALSE, segment.color = "grey60") +
      stat_ellipse(aes(group = Group), type = "norm",
                   linetype = "dashed", alpha = 0.6, linewidth = 0.7) +
      scale_color_manual(values = group_palette) +
      annotate("label", x = Inf, y = -Inf,
               label      = annot,
               hjust      = 1.03, vjust = -0.15,
               size       = 2.8, fill = "white", alpha = 0.88,
               label.size = 0.35, family = "mono") +
      labs(
        title    = sprintf("PCoA - %s - %s", ds_name, dist_name),
        subtitle = sprintf("Rarefied to %d counts | 999 permutations", rare_depth),
        x = sprintf("PCoA Axis 1 (%.1f%%)", var_pcoa[1]),
        y = sprintf("PCoA Axis 2 (%.1f%%)", var_pcoa[2]),
        color = group_col_name
      ) +
      theme_bw(base_size = 13) +
      theme(
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey40", size = 10),
        legend.position = "right"
      )
    
    pcoa_plots[[dist_name]] <- p_pcoa
    print(p_pcoa)
  }
  
  # PCoA panel
  p_pcoa_panel <- wrap_plots(pcoa_plots, ncol = 2) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = paste("PCoA - Beta Diversity -", ds_name),
      theme = theme(plot.title = element_text(face = "bold", size = 15))
    ) &
    theme(legend.position = "bottom")
  
  print(p_pcoa_panel)
  
  # 9.6 NMDS
  nmds_plots  <- list()
  nmds_stress <- list()
  
  for (dist_name in names(dist_list)) {
    
    set.seed(42)
    nmds_res <- tryCatch({
      metaMDS(dist_list[[dist_name]],
              k = 2, trymax = 100,
              trace = FALSE, autotransform = FALSE)
    }, error = function(e) {
      cat(paste("    AVISO: NMDS falhou para", dist_name, ":", e$message, "\n"))
      NULL
    })
    
    if (is.null(nmds_res)) next
    
    stress_val <- nmds_res$stress
    nmds_stress[[dist_name]] <- stress_val
    
    stress_label <- ifelse(stress_val < 0.10, "excelente",
                           ifelse(stress_val < 0.20, "bom",
                                  ifelse(stress_val < 0.30, "aceitavel", "ruim")))
    
    cat(sprintf("    NMDS [%s]: stress=%.4f (%s)\n",
                dist_name, stress_val, stress_label))
    
    nmds_scores <- as.data.frame(scores(nmds_res, display = "sites"))
    nmds_scores$Sample <- rownames(nmds_scores)
    nmds_scores$Group  <- meta_beta$Group[match(nmds_scores$Sample, meta_beta$Sample)]
    
    p_nmds <- ggplot(nmds_scores,
                     aes(x = NMDS1, y = NMDS2,
                         color = Group, label = Sample)) +
      geom_point(size = 3.5, alpha = 0.85) +
      geom_text_repel(size = 2.5, max.overlaps = 25,
                      show.legend = FALSE, segment.color = "grey60") +
      stat_ellipse(aes(group = Group), type = "norm",
                   linetype = "dashed", alpha = 0.6, linewidth = 0.7) +
      scale_color_manual(values = group_palette) +
      annotate("text", x = -Inf, y = Inf,
               label    = sprintf("Stress: %.4f (%s)", stress_val, stress_label),
               hjust    = -0.1, vjust = 1.4,
               size     = 3.5, color = "grey30", fontface = "italic") +
      labs(
        title    = sprintf("NMDS - %s - %s", ds_name, dist_name),
        subtitle = sprintf("Stress=%.4f (%s) | Rarefied to %d",
                           stress_val, stress_label, rare_depth),
        x = "NMDS1", y = "NMDS2", color = group_col_name
      ) +
      theme_bw(base_size = 13) +
      theme(
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey40", size = 10),
        legend.position = "right"
      )
    
    nmds_plots[[dist_name]] <- p_nmds
    print(p_nmds)
  }
  
  # NMDS panel
  p_nmds_panel <- NULL
  if (length(nmds_plots) > 0) {
    p_nmds_panel <- wrap_plots(nmds_plots, ncol = 2) +
      plot_layout(guides = "collect") +
      plot_annotation(
        title = paste("NMDS - Beta Diversity -", ds_name),
        theme = theme(plot.title = element_text(face = "bold", size = 15))
      ) &
      theme(legend.position = "bottom")
    
    print(p_nmds_panel)
  }
  
  # 9.7 Heatmap — Top features mais variáveis
  n_top_heat <- min(50, ncol(rel_mat))
  
  feat_vars <- apply(rel_mat, 2, var)
  top_feats <- names(sort(feat_vars, decreasing = TRUE))[1:n_top_heat]
  
  heat_mat   <- t(rel_mat[, top_feats])
  heat_mat_z <- t(scale(t(heat_mat)))
  heat_mat_z <- pmax(pmin(heat_mat_z, 3), -3)
  heat_mat_z[is.nan(heat_mat_z)] <- 0
  
  # Clusterização
  hclust_cols <- tryCatch(
    hclust(dist_bray, method = "ward.D2"),
    error = function(e) hclust(dist(t(heat_mat_z)), method = "ward.D2")
  )
  hclust_rows <- hclust(dist(heat_mat_z), method = "ward.D2")
  
  col_fun <- colorRamp2(
    c(-3, -1.5, 0, 1.5, 3),
    c("#2166AC", "#92C5DE", "white", "#F4A582", "#D6604D")
  )
  
  col_anno <- HeatmapAnnotation(
    Group = meta_beta$Group,
    col   = list(Group = group_palette),
    annotation_name_gp = gpar(fontsize = 11, fontface = "bold")
  )
  
  row_fontsize <- ifelse(n_top_heat > 40, 5.5, ifelse(n_top_heat > 25, 7, 8.5))
  
  ht <- Heatmap(
    heat_mat_z,
    name              = "Z-score",
    col               = col_fun,
    top_annotation    = col_anno,
    cluster_rows      = hclust_rows,
    cluster_columns   = hclust_cols,
    show_column_names = TRUE,
    column_names_gp   = gpar(fontsize = 7),
    show_row_names    = TRUE,
    row_names_gp      = gpar(fontsize = row_fontsize),
    row_names_side    = "left",
    row_dend_side     = "right",
    column_title      = sprintf(
      "%s - Top %d Variable Features (Z-score, Bray-Curtis clustering)",
      ds_name, n_top_heat
    ),
    column_title_gp   = gpar(fontsize = 12, fontface = "bold"),
    heatmap_legend_param = list(
      title    = "Z-score",
      at       = c(-3, -1.5, 0, 1.5, 3),
      title_gp = gpar(fontface = "bold")
    ),
    border = TRUE
  )
  
  draw(ht, merge_legend = TRUE)
  
  # 9.8 Armazenar resultados
  permanova_tbl <- do.call(rbind, lapply(names(stats_by_dist), function(dn) {
    s <- stats_by_dist[[dn]]
    data.frame(
      Dataset  = ds_name,
      Distance = dn,
      R2       = s$permanova_R2,
      F_stat   = s$permanova_F,
      p_value  = s$permanova_p,
      Sig      = s$permanova_sig,
      stringsAsFactors = FALSE
    )
  }))
  
  anosim_tbl <- do.call(rbind, lapply(names(stats_by_dist), function(dn) {
    s <- stats_by_dist[[dn]]
    data.frame(
      Dataset  = ds_name,
      Distance = dn,
      R_stat   = s$anosim_R,
      p_value  = s$anosim_p,
      Strength = ifelse(s$anosim_R > 0.75, "strong",
                        ifelse(s$anosim_R > 0.50, "moderate",
                               ifelse(s$anosim_R > 0.25, "weak", "negligible"))),
      Sig      = s$anosim_sig,
      stringsAsFactors = FALSE
    )
  }))
  
  nmds_stress_tbl <- if (length(nmds_stress) > 0) {
    data.frame(
      Dataset  = ds_name,
      Distance = names(nmds_stress),
      Stress   = round(unlist(nmds_stress), 4),
      Quality  = sapply(unlist(nmds_stress), function(s) {
        ifelse(s < 0.10, "Excellent",
               ifelse(s < 0.20, "Good",
                      ifelse(s < 0.30, "Acceptable", "Poor")))
      }),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(Dataset = character(), Distance = character(),
               Stress = numeric(), Quality = character(),
               stringsAsFactors = FALSE)
  }
  
  beta_stats_all[[ds_name]] <- list(
    permanova   = permanova_tbl,
    anosim      = anosim_tbl,
    betadisper  = betadisp_results,
    nmds_stress = nmds_stress_tbl
  )
  
  beta_plots_all[[ds_name]] <- list(
    pcoa_plots = pcoa_plots,
    pcoa_panel = p_pcoa_panel,
    nmds_plots = nmds_plots,
    nmds_panel = p_nmds_panel,
    heatmap    = ht
  )
  
  cat(paste("\n  OK:", ds_name, "beta diversity completa!\n\n"))
}

# 9.9 Resumo consolidado Beta Diversity
cat("\n", rep("=", 70), "\n")
cat("=== RESUMO CONSOLIDADO - BETA DIVERSITY ===\n")
cat(rep("=", 70), "\n\n")

all_permanova <- do.call(rbind, lapply(beta_stats_all, function(x) x$permanova))
all_anosim    <- do.call(rbind, lapply(beta_stats_all, function(x) x$anosim))
all_betadisp  <- do.call(rbind, lapply(beta_stats_all, function(x) x$betadisper))
all_nmds_str  <- do.call(rbind, lapply(beta_stats_all, function(x) x$nmds_stress))

cat("PERMANOVA:\n")
print(all_permanova)
cat("\nANOSIM:\n")
print(all_anosim)
cat("\nBetadisper:\n")
print(all_betadisp)
cat("\nNMDS Stress:\n")
print(all_nmds_str)

# ═══════════════════════════════════════════════════════════════════════════════
# ███████████████████████████████████████████████████████████████████████████████
# FASE 10 — EXPORT: PLOTS + TABELAS EXCEL
# ███████████████████████████████████████████████████████████████████████████████
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n", rep("=", 70), "\n")
cat("=== FASE 10: EXPORTACAO DE RESULTADOS ===\n")
cat(rep("=", 70), "\n\n")

# 10.1 Criar estrutura de diretórios
for (ds_name in names(datasets_grouped_list)) {
  dirs_to_create <- c(
    file.path(output_root, ds_name, "plots", "pca"),
    file.path(output_root, ds_name, "plots", "alpha_diversity"),
    file.path(output_root, ds_name, "plots", "beta_diversity"),
    file.path(output_root, ds_name, "tables")
  )
  for (d in dirs_to_create) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
dir.create(file.path(output_root, "consolidated"), recursive = TRUE, showWarnings = FALSE)

# 10.2 Helper para formatar sheet Excel
add_fmt_sheet <- function(wb, sheet_name, data, hdr_color = "#4472C4") {
  addWorksheet(wb, sheetName = sheet_name)
  
  hdr_style <- createStyle(
    fontColour = "#FFFFFF", fgFill = hdr_color,
    fontName = "Calibri", fontSize = 11,
    textDecoration = "bold", halign = "center", valign = "center",
    border = "TopBottomLeftRight", borderColour = "#FFFFFF", wrapText = TRUE
  )
  body_style <- createStyle(
    fontName = "Calibri", fontSize = 10,
    border = "TopBottomLeftRight", borderColour = "#BFBFBF"
  )
  alt_style <- createStyle(
    fontName = "Calibri", fontSize = 10, fgFill = "#F2F7FF",
    border = "TopBottomLeftRight", borderColour = "#BFBFBF"
  )
  
  writeData(wb, sheet = sheet_name, x = data, headerStyle = hdr_style)
  
  nr <- nrow(data)
  nc <- ncol(data)
  if (nr > 0) {
    addStyle(wb, sheet_name, style = body_style,
             rows = 2:(nr + 1), cols = 1:nc, gridExpand = TRUE)
    even_rows <- seq(3, nr + 1, by = 2)
    if (length(even_rows) > 0) {
      addStyle(wb, sheet_name, style = alt_style,
               rows = even_rows, cols = 1:nc, gridExpand = TRUE)
    }
  }
  
  setColWidths(wb, sheet_name, cols = 1:nc, widths = "auto")
  freezePane(wb, sheet_name, firstRow = TRUE)
  invisible(wb)
}

# 10.3 Exportar plots por dataset
cat("Exportando plots...\n\n")

for (ds_name in names(datasets_grouped_list)) {
  
  cat(paste("  [", ds_name, "]\n"))
  
  # PCA plots
  d_pca <- file.path(output_root, ds_name, "plots", "pca")
  if (!is.null(pca_results_grouped[[ds_name]])) {
    pca_res <- pca_results_grouped[[ds_name]]
    save_gg(pca_res$plot_pc12, file.path(d_pca, "pca_PC1_PC2"), width = 9, height = 7)
    if (!is.null(pca_res$plot_pc13)) {
      save_gg(pca_res$plot_pc13, file.path(d_pca, "pca_PC1_PC3"), width = 9, height = 7)
    }
  }
  
  # Alpha plots
  d_alpha <- file.path(output_root, ds_name, "plots", "alpha_diversity")
  if (!is.null(alpha_plots_all[[ds_name]])) {
    ap <- alpha_plots_all[[ds_name]]
    save_gg(ap$depth,    file.path(d_alpha, "sequencing_depth"),     width = 13, height = 6)
    save_gg(ap$facet,    file.path(d_alpha, "alpha_all_metrics"),    width = 12, height = 10)
    save_gg(ap$panel,    file.path(d_alpha, "alpha_2x2_panel"),      width = 14, height = 12)
    save_gg(ap$observed, file.path(d_alpha, "alpha_observed"),       width = 7,  height = 6)
    save_gg(ap$chao1,    file.path(d_alpha, "alpha_chao1"),          width = 7,  height = 6)
    save_gg(ap$shannon,  file.path(d_alpha, "alpha_shannon"),        width = 7,  height = 6)
    save_gg(ap$simpson,  file.path(d_alpha, "alpha_simpson"),        width = 7,  height = 6)
  }
  
  # Beta plots
  d_beta <- file.path(output_root, ds_name, "plots", "beta_diversity")
  if (!is.null(beta_plots_all[[ds_name]])) {
    bp <- beta_plots_all[[ds_name]]
    
    # PCoA
    for (dist_name in names(bp$pcoa_plots)) {
      fn <- gsub("[^a-zA-Z0-9]", "_", dist_name)
      save_gg(bp$pcoa_plots[[dist_name]],
              file.path(d_beta, paste0("pcoa_", fn)), width = 9, height = 7)
    }
    save_gg(bp$pcoa_panel, file.path(d_beta, "pcoa_panel"), width = 16, height = 7)
    
    # NMDS
    for (dist_name in names(bp$nmds_plots)) {
      fn <- gsub("[^a-zA-Z0-9]", "_", dist_name)
      save_gg(bp$nmds_plots[[dist_name]],
              file.path(d_beta, paste0("nmds_", fn)), width = 9, height = 7)
    }
    if (!is.null(bp$nmds_panel)) {
      save_gg(bp$nmds_panel, file.path(d_beta, "nmds_panel"), width = 16, height = 7)
    }
    
    # Heatmap
    tryCatch({
      png(file.path(d_beta, "heatmap_top_features.png"),
          width = 14, height = 12, units = "in", res = 300, bg = "white")
      draw(bp$heatmap, merge_legend = TRUE)
      dev.off()
      cat(paste("    Salvo: heatmap_top_features.png\n"))
    }, error = function(e) {
      cat(paste("    AVISO: Heatmap export falhou:", e$message, "\n"))
    })
  }
}

# 10.4 Exportar tabelas Excel por dataset
cat("\nExportando tabelas Excel...\n\n")

for (ds_name in names(datasets_grouped_list)) {
  
  tbl_dir <- file.path(output_root, ds_name, "tables")
  
  # Alpha Excel
  wb_alpha <- createWorkbook()
  add_fmt_sheet(wb_alpha, "Descriptive",
                alpha_stats_all[[ds_name]]$descriptive, hdr_color = "#2E75B6")
  add_fmt_sheet(wb_alpha, "Kruskal_Wallis",
                alpha_stats_all[[ds_name]]$kruskal_wallis, hdr_color = "#2E75B6")
  add_fmt_sheet(wb_alpha, "Wilcoxon",
                alpha_stats_all[[ds_name]]$wilcoxon, hdr_color = "#2E75B6")
  add_fmt_sheet(wb_alpha, "Per_Sample",
                alpha_results_all[[ds_name]], hdr_color = "#2E75B6")
  saveWorkbook(wb_alpha,
               file.path(tbl_dir, paste0(ds_name, "_alpha_diversity.xlsx")),
               overwrite = TRUE)
  cat(paste("  ", ds_name, "_alpha_diversity.xlsx\n"))
  
  # Beta Excel
  wb_beta <- createWorkbook()
  add_fmt_sheet(wb_beta, "PERMANOVA",
                beta_stats_all[[ds_name]]$permanova, hdr_color = "#375623")
  add_fmt_sheet(wb_beta, "ANOSIM",
                beta_stats_all[[ds_name]]$anosim, hdr_color = "#375623")
  add_fmt_sheet(wb_beta, "Betadisper",
                beta_stats_all[[ds_name]]$betadisper, hdr_color = "#375623")
  add_fmt_sheet(wb_beta, "NMDS_Stress",
                beta_stats_all[[ds_name]]$nmds_stress, hdr_color = "#375623")
  saveWorkbook(wb_beta,
               file.path(tbl_dir, paste0(ds_name, "_beta_diversity.xlsx")),
               overwrite = TRUE)
  cat(paste("  ", ds_name, "_beta_diversity.xlsx\n"))
}

# 10.5 Tabela MASTER consolidada
wb_master <- createWorkbook()
add_fmt_sheet(wb_master, "Alpha_Descriptive",    all_desc,       hdr_color = "#2E75B6")
add_fmt_sheet(wb_master, "Alpha_Kruskal_Wallis", all_kw,         hdr_color = "#2E75B6")
add_fmt_sheet(wb_master, "Alpha_Wilcoxon",       all_wx,         hdr_color = "#2E75B6")
add_fmt_sheet(wb_master, "Beta_PERMANOVA",       all_permanova,  hdr_color = "#375623")
add_fmt_sheet(wb_master, "Beta_ANOSIM",          all_anosim,     hdr_color = "#375623")
add_fmt_sheet(wb_master, "Beta_Betadisper",      all_betadisp,   hdr_color = "#375623")
add_fmt_sheet(wb_master, "Beta_NMDS_Stress",     all_nmds_str,   hdr_color = "#375623")

# PCA variance summary
if (exists("variance_summary") && !is.null(variance_summary)) {
  add_fmt_sheet(wb_master, "PCA_Variance", variance_summary, hdr_color = "#7030A0")
}

# Rarefação info
rare_info <- do.call(rbind, rarefaction_info_all)
add_fmt_sheet(wb_master, "Rarefaction_Info", rare_info, hdr_color = "#7030A0")

saveWorkbook(wb_master,
             file.path(output_root, "consolidated", "MASTER_diversity_all_datasets.xlsx"),
             overwrite = TRUE)
cat("\n  MASTER_diversity_all_datasets.xlsx\n")

# 10.6 Session info
session_path <- file.path(output_root, "session_info.txt")
sink(session_path)
cat("=== FUNCTIONAL DIVERSITY ANALYSIS — SESSION INFO ===\n")
cat(sprintf("Generated: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Data dir:  %s\n", data_dir))
cat(sprintf("Groups:    %s\n", paste(selected_groups, collapse = ", ")))
cat(sprintf("Ages:      %s\n", ifelse(is.null(selected_ages), "ALL",
                                      paste(selected_ages, collapse = ", "))))
cat("\n")
print(sessionInfo())
sink()

# 10.7 Sumário final
cat("\n", rep("=", 70), "\n")
cat("  PIPELINE COMPLETO - EXPORT FINALIZADO\n")
cat(rep("=", 70), "\n\n")
cat(sprintf("  Root: %s/\n\n", output_root))
cat("  Estrutura por dataset:\n")

for (ds_name in names(datasets_grouped_list)) {
  cat(sprintf("  %s/\n", ds_name))
  cat("    plots/\n")
  cat("      pca/                 PCA PC1xPC2, PC1xPC3\n")
  cat("      alpha_diversity/     depth, facet, 2x2, individual metrics\n")
  cat("      beta_diversity/      PCoA, NMDS (Bray-Curtis + Jaccard), heatmap\n")
  cat("    tables/\n")
  cat(sprintf("      %s_alpha_diversity.xlsx\n", ds_name))
  cat(sprintf("      %s_beta_diversity.xlsx\n", ds_name))
}

cat("\n  consolidated/\n")
cat("    MASTER_diversity_all_datasets.xlsx   <- tudo em um arquivo\n")
cat("  session_info.txt\n\n")
cat("  Formato: PNG 300dpi | Excel formatado com headers coloridos\n")
cat(rep("=", 70), "\n")
cat("PIPELINE FASES 1-10 CONCLUIDO COM SUCESSO!\n")
cat(rep("=", 70), "\n")
