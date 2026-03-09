#Step -by step
# ============================================================
# MICROBIOME ANALYSIS — COMPLETE REVISED CODE v001
# ============================================================

# ─────────────────────────────────────────────────────────────────────────────
# PART 0 — DIR
# ─────────────────────────────────────────────────────────────────────────────
setwd("YOUR_DIR")
# ─────────────────────────────────────────────────────────────────────────────
# PART 0 — INSTALLATION
# ─────────────────────────────────────────────────────────────────────────────

# Verifica e instala 'BiocManager' se necessário
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Lista de pacotes Bioconductor
bioc_packages <- c(
  "phyloseq", "microbiome", "DESeq2", "ALDEx2",
  "lefser", "KEGGREST", "clusterProfiler",
  "MicrobiomeStat", "Maaslin2", "ComplexHeatmap"
)

# Instala pacotes Bioconductor se não estiverem presentes
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# Lista de pacotes CRAN
cran_packages <- c(
  "tidyverse", "vegan", "ggplot2", "ggrepel",
  "reshape2", "RColorBrewer", "pheatmap",
  "ggpubr", "scales", "patchwork", "ape", "circlize"
)

# Instala pacotes CRAN se não estiverem presentes
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
install.packages("openxlsx")
# ─────────────────────────────────────────────────────────────────────────────
# PART 1 — LIBRARIES
# ─────────────────────────────────────────────────────────────────────────────

# Carrega as bibliotecas
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(microbiome)
library(Maaslin2)
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(ape)
library(openxlsx)
# ─────────────────────────────────────────────────────────────────────────────
# PART 2 — DATA READING AND ID STANDARDIZATION
# ─────────────────────────────────────────────────────────────────────────────

fix_id   <- function(x) sub("^X(?=\\d)", "", x, perl = TRUE)
tax_cols <- c("k", "p", "c", "o", "f", "g", "s")

# 2.1  Reading
tax_counts_df <- read_csv(
  "ANH/Simphyome_Analysis/Mick_example_Files/hav2101_all_samples_counts_wide_tsv_converted_to_csv.csv"
)
metadata_df <- read_csv(
  "ANH/Simphyome_Analysis/Mick_example_Files/HAV2101_sample_metadata.csv"
)
ko_counts_df <- read.table(
  "ANH/Simphyome_Analysis/Mick_example_Files/hav2101_ko_count_table.txt",
  header = TRUE, sep = "\t", dec = "."
)

# 2.2  Standardize column names
colnames(tax_counts_df) <- fix_id(colnames(tax_counts_df))
colnames(ko_counts_df)  <- str_remove_all(colnames(ko_counts_df), "^X|\\.counts\\.$")

# 2.3  Standardize metadata IDs
metadata_df <- metadata_df %>%
  mutate(Sample_ID_FM_Pipeline = fix_id(as.character(Sample_ID_FM_Pipeline)))

# ─────────────────────────────────────────────────────────────────────────────
# PART 3 — INITIAL FILTERING
# ─────────────────────────────────────────────────────────────────────────────

metadata_filtered <- metadata_df %>%
  filter(Group %in% c("Cp Challenge Control", "Symphiome"))

ids_meta <- metadata_filtered$Sample_ID_FM_Pipeline
ids_tax  <- setdiff(colnames(tax_counts_df), tax_cols)
ids_ko   <- setdiff(colnames(ko_counts_df),  "ko")

cat("\nSamples in KO but missing in metadata:\n")
print(setdiff(ids_ko,  ids_meta))
cat("\nSamples in TAX but missing in metadata:\n")
print(setdiff(ids_tax, ids_meta))

ids_common <- Reduce(intersect, list(ids_meta, ids_tax, ids_ko))
cat(sprintf("\nCommon samples: %d\n", length(ids_common)))

tax_counts_df_filtrado <- tax_counts_df %>%
  select(all_of(tax_cols), all_of(ids_common))

ko_counts_df_filtrado <- ko_counts_df %>%
  select(ko, all_of(ids_common))

metadata_filtered <- metadata_filtered %>%
  filter(Sample_ID_FM_Pipeline %in% ids_common) %>%
  arrange(match(Sample_ID_FM_Pipeline, ids_common))
# ─────────────────────────────────────────────────────────────────────────────
# PART 4 — QC: HOST CONTAMINATION
# ─────────────────────────────────────────────────────────────────────────────

# 4.1  Host reads
host_counts <- tax_counts_df_filtrado %>%
  filter(s == "Gallus_gallus") %>%
  select(-all_of(tax_cols)) %>%
  pivot_longer(everything(),
               names_to  = "Sample_ID_FM_Pipeline",
               values_to = "Host_Reads")

# 4.2  Total reads per sample
total_reads <- tax_counts_df_filtrado %>%
  select(-all_of(tax_cols)) %>%
  summarise(across(everything(), sum)) %>%
  pivot_longer(everything(),
               names_to  = "Sample_ID_FM_Pipeline",
               values_to = "Total_Reads")

# 4.3  Contamination percentage
host_contamination_df <- inner_join(host_counts, total_reads,
                                    by = "Sample_ID_FM_Pipeline") %>%
  filter(Total_Reads > 0) %>%
  mutate(Host_Pct = (Host_Reads / Total_Reads) * 100)

# 4.4  Plot — all samples
p_host_all <- ggplot(
  host_contamination_df,
  aes(x = reorder(Sample_ID_FM_Pipeline, -Host_Pct), y = Host_Pct)
) +
  geom_segment(aes(xend = Sample_ID_FM_Pipeline, yend = 0), color = "grey70") +
  geom_point(color = "darkred", size = 3) +
  labs(title = "Host Contamination - All Samples",
       x = "Sample", y = "Contamination (%)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

print(p_host_all)

# 4.5  Plot — samples > 0% only
p_host_nonzero <- host_contamination_df %>%
  filter(Host_Pct > 0) %>%
  ggplot(aes(x = reorder(Sample_ID_FM_Pipeline, -Host_Pct), y = Host_Pct)) +
  geom_segment(aes(xend = Sample_ID_FM_Pipeline, yend = 0), color = "grey70") +
  geom_point(color = "darkred", size = 3) +
  labs(title = "Host Contamination (> 0%)",
       x = "Sample", y = "Contamination (%)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

print(p_host_nonzero)

# 4.6  Filter samples and remove host row
host_threshold <- 20

samples_to_keep <- host_contamination_df %>%
  filter(Host_Pct < host_threshold) %>%
  pull(Sample_ID_FM_Pipeline)

tax_counts_final_qc <- tax_counts_df_filtrado %>%
  filter(s != "Gallus_gallus") %>%
  select(all_of(tax_cols), all_of(samples_to_keep))

ko_counts_final_qc <- ko_counts_df_filtrado %>%
  select(ko, all_of(samples_to_keep))

metadata_final_qc <- metadata_filtered %>%
  filter(Sample_ID_FM_Pipeline %in% samples_to_keep)

cat(sprintf(
  "\nThreshold: %d%% | Removed: %d | Remaining: %d\n",
  host_threshold,
  length(ids_common) - length(samples_to_keep),
  length(samples_to_keep)
))

# ─────────────────────────────────────────────────────────────────────────────
# PART 4.7 — FILTER LOW-ABUNDANCE / LOW-PREVALENCE FEATURES
#            Keep species/KOs with ≥ 0.01% rel. abundance in ≥ 10% of samples
# ─────────────────────────────────────────────────────────────────────────────

min_rel_abund  <- 0.0001   # 0.01% relative abundance
min_prevalence <- 0.10     # present in at least 10% of samples

# ── Taxonomic filtering ───────────────────────────────────────────────────────

# Step 1: extract count matrix only (no taxonomy columns)
tax_count_mat <- tax_counts_final_qc %>%
  select(-all_of(tax_cols)) %>%
  as.matrix()

# Step 2: per-sample relative abundance (column-wise proportions)
tax_rel_mat <- sweep(tax_count_mat, 2, colSums(tax_count_mat), "/")

# Step 3: for each species (row), fraction of samples exceeding threshold
tax_prevalence <- rowMeans(tax_rel_mat >= min_rel_abund, na.rm = TRUE)

# Step 4: flag which rows pass
tax_rows_pass <- tax_prevalence >= min_prevalence

# Step 5: apply filter, keeping taxonomy columns
tax_counts_final_qc <- tax_counts_final_qc[tax_rows_pass, ]

cat(sprintf(
  "\n[TAX] Species before prevalence filter: %d | After: %d | Removed: %d\n",
  length(tax_rows_pass),
  sum(tax_rows_pass),
  sum(!tax_rows_pass)
))

# ── KO filtering ──────────────────────────────────────────────────────────────

# Step 1: extract count matrix only (no ko column)
ko_count_mat <- ko_counts_final_qc %>%
  select(-ko) %>%
  as.matrix()

# Step 2: per-sample relative abundance
ko_rel_mat <- sweep(ko_count_mat, 2, colSums(ko_count_mat), "/")

# Step 3: fraction of samples exceeding threshold per KO
ko_prevalence <- rowMeans(ko_rel_mat >= min_rel_abund, na.rm = TRUE)

# Step 4: flag passing rows
ko_rows_pass <- ko_prevalence >= min_prevalence

# Step 5: apply filter, keeping ko column
ko_counts_final_qc <- ko_counts_final_qc[ko_rows_pass, ]

cat(sprintf(
  "[KO]  Features before prevalence filter: %d | After: %d | Removed: %d\n",
  length(ko_rows_pass),
  sum(ko_rows_pass),
  sum(!ko_rows_pass)
))

# ── Summary ───────────────────────────────────────────────────────────────────

cat(sprintf(
  "\nFilter thresholds: ≥ %.2f%% relative abundance in ≥ %.0f%% of samples\n",
  min_rel_abund * 100,
  min_prevalence * 100
))

#-----------------------------------------------------------------------
#PCA
# ─────────────────────────────────────────────────────────────────────────────
# PART 5 — PCA: SAMPLE SWAP CHECK
# ─────────────────────────────────────────────────────────────────────────────

library(ggrepel)   # for non-overlapping sample labels

# 5.1  Prepare species-level count matrix ─────────────────────────────────────

pca_matrix <- tax_counts_final_qc %>%
  filter(!is.na(s)) %>%                      # keep only rows with species annotation
  select(-all_of(tax_cols)) %>%              # drop taxonomy columns
  filter(rowSums(.) > 0)                     # remove all-zero rows

# Transpose → samples as rows, species as columns
pca_matrix_t <- t(pca_matrix)

# 5.2  CLR transformation (appropriate for compositional data) ────────────────

pca_clr <- pca_matrix_t + 0.5                                     # pseudocount
pca_clr <- sweep(pca_clr, 1, rowSums(pca_clr), "/")              # relative abundances
pca_clr <- log(pca_clr) - rowMeans(log(pca_clr))                 # CLR

# 5.3  Run PCA ─────────────────────────────────────────────────────────────────

pca_result   <- prcomp(pca_clr, center = TRUE, scale. = FALSE)
var_exp      <- summary(pca_result)$importance[2, ] * 100        # % variance explained

# 5.4  Build plotting data frame ───────────────────────────────────────────────

pca_df <- as.data.frame(pca_result$x[, 1:3]) %>%
  rownames_to_column("Sample_ID_FM_Pipeline") %>%
  left_join(
    metadata_final_qc %>% select(Sample_ID_FM_Pipeline, Group),
    by = "Sample_ID_FM_Pipeline"
  )

# 5.5  PC1 vs PC2 ──────────────────────────────────────────────────────────────

p_pca_12 <- ggplot(pca_df,
                   aes(x = PC1, y = PC2,
                       color = Group,
                       label = Sample_ID_FM_Pipeline)) +
  geom_point(size = 3.5, alpha = 0.85) +
  geom_text_repel(
    size         = 3,
    max.overlaps = 30,
    show.legend  = FALSE,
    segment.color = "grey60"
  ) +
  stat_ellipse(
    aes(group = Group),
    type     = "norm",
    linetype = "dashed",
    alpha    = 0.5,
    linewidth = 0.7
  ) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title    = "PCA — Sample Swap Check (PC1 vs PC2)",
    subtitle = "Species-level · CLR transformation",
    x        = sprintf("PC1 (%.1f%%)", var_exp[1]),
    y        = sprintf("PC2 (%.1f%%)", var_exp[2]),
    color    = "Group"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40"),
    legend.position = "right"
  )

print(p_pca_12)

# 5.6  PC1 vs PC3 (optional — extra swap check) ────────────────────────────────

p_pca_13 <- ggplot(pca_df,
                   aes(x = PC1, y = PC3,
                       color = Group,
                       label = Sample_ID_FM_Pipeline)) +
  geom_point(size = 3.5, alpha = 0.85) +
  geom_text_repel(
    size         = 3,
    max.overlaps = 30,
    show.legend  = FALSE,
    segment.color = "grey60"
  ) +
  stat_ellipse(
    aes(group = Group),
    type     = "norm",
    linetype = "dashed",
    alpha    = 0.5,
    linewidth = 0.7
  ) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title    = "PCA — Sample Swap Check (PC1 vs PC3)",
    subtitle = "Species-level · CLR transformation",
    x        = sprintf("PC1 (%.1f%%)", var_exp[1]),
    y        = sprintf("PC3 (%.1f%%)", var_exp[3]),
    color    = "Group"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40"),
    legend.position = "right"
  )

print(p_pca_13)

# 5.7  Scree plot (optional — how many PCs matter?) ────────────────────────────

scree_df <- data.frame(
  PC      = seq_along(var_exp),
  Var_Pct = var_exp
)

p_scree <- ggplot(scree_df[1:15, ], aes(x = PC, y = Var_Pct)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_line(color = "darkred", linewidth = 0.8) +
  geom_point(color = "darkred", size = 2) +
  labs(
    title = "Scree Plot — Variance Explained per PC",
    x     = "Principal Component",
    y     = "Variance Explained (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

print(p_scree)

# ─────────────────────────────────────────────────────────────────────────────
# 5.8  Biplot — samples + species loadings
# ─────────────────────────────────────────────────────────────────────────────

n_arrows <- 15   # number of top-contributing species to display as arrows

# Retrieve species names matching pca_matrix row order
# NOTE: must apply identical filters as in 5.1 to guarantee row alignment
species_labels_biplot <- tax_counts_final_qc %>%
  filter(!is.na(s)) %>%
  mutate(count_sum = rowSums(select(., -all_of(tax_cols)))) %>%
  filter(count_sum > 0) %>%
  mutate(species_label = make.unique(paste(g, s, sep = "|"))) %>%
  pull(species_label)

# Assign species names as rownames of the loading matrix
rownames(pca_result$rotation) <- species_labels_biplot

# ─────────────────────────────────────────────────────────────────────────────
# 5.8.1  Helper — build scaled loadings data frame for any pair of PCs
# ─────────────────────────────────────────────────────────────────────────────

get_scaled_loadings <- function(rotation, scores_df,
                                pc_x, pc_y,
                                n_top, scale_pct = 0.65) {
  
  load_df <- as.data.frame(rotation[, c(pc_x, pc_y)]) %>%
    setNames(c("PCx", "PCy")) %>%
    rownames_to_column("species") %>%
    mutate(
      species_short = sub("^[^|]*\\|", "", species),   # keep part after "|"
      loading_mag   = sqrt(PCx^2 + PCy^2)
    ) %>%
    arrange(desc(loading_mag)) %>%
    head(n_top)
  
  # Scale so arrows reach scale_pct of the sample-score range
  score_range   <- max(abs(scores_df[, c(pc_x, pc_y)]))
  loading_range <- max(abs(load_df[, c("PCx", "PCy")]))
  sf            <- score_range / loading_range * scale_pct
  
  load_df %>%
    mutate(PCx_sc = PCx * sf,
           PCy_sc = PCy * sf)
}

# ─────────────────────────────────────────────────────────────────────────────
# 5.8.2  Helper — build biplot ggplot
# ─────────────────────────────────────────────────────────────────────────────

make_biplot <- function(scores_df, loadings_scaled,
                        pc_x_col, pc_y_col,
                        x_label, y_label,
                        title, subtitle) {
  
  ggplot() +
    # ── Sample points ──────────────────────────────────────────────────────
    geom_point(
      data  = scores_df,
      aes(x = .data[[pc_x_col]], y = .data[[pc_y_col]], color = Group),
      size  = 3.2, alpha = 0.82
    ) +
    geom_text_repel(
      data          = scores_df,
      aes(x = .data[[pc_x_col]], y = .data[[pc_y_col]],
          color = Group, label = Sample_ID_FM_Pipeline),
      size          = 2.6,
      max.overlaps  = 20,
      show.legend   = FALSE,
      segment.color = "grey60",
      segment.size  = 0.3
    ) +
    # ── Confidence ellipses ────────────────────────────────────────────────
    stat_ellipse(
      data      = scores_df,
      aes(x = .data[[pc_x_col]], y = .data[[pc_y_col]],
          group = Group, color = Group),
      type      = "norm",
      linetype  = "dashed",
      alpha     = 0.5,
      linewidth = 0.7
    ) +
    # ── Loading arrows ─────────────────────────────────────────────────────
    geom_segment(
      data      = loadings_scaled,
      aes(x = 0, y = 0, xend = PCx_sc, yend = PCy_sc),
      arrow     = arrow(length = unit(0.22, "cm"), type = "closed"),
      color     = "grey20",
      linewidth = 0.45,
      alpha     = 0.80
    ) +
    # ── Arrow labels ───────────────────────────────────────────────────────
    geom_text_repel(
      data          = loadings_scaled,
      aes(x = PCx_sc, y = PCy_sc, label = species_short),
      size          = 2.5,
      color         = "grey10",
      fontface      = "italic",
      max.overlaps  = 40,
      segment.color = "grey65",
      segment.size  = 0.25,
      box.padding   = 0.35,
      force         = 1.5
    ) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title    = title,
      subtitle = subtitle,
      x        = x_label,
      y        = y_label,
      color    = "Group"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title      = element_text(face = "bold"),
      plot.subtitle   = element_text(color = "grey40", size = 10),
      legend.position = "right"
    )
}

# ─────────────────────────────────────────────────────────────────────────────
# 5.8.3  Biplot — PC1 vs PC2
# ─────────────────────────────────────────────────────────────────────────────

load_scaled_12 <- get_scaled_loadings(
  rotation  = pca_result$rotation,
  scores_df = pca_df,
  pc_x      = "PC1",
  pc_y      = "PC2",
  n_top     = n_arrows
)

p_biplot_12 <- make_biplot(
  scores_df       = pca_df,
  loadings_scaled = load_scaled_12,
  pc_x_col        = "PC1",
  pc_y_col        = "PC2",
  x_label         = sprintf("PC1 (%.1f%%)", var_exp[1]),
  y_label         = sprintf("PC2 (%.1f%%)", var_exp[2]),
  title           = "PCA Biplot — PC1 vs PC2",
  subtitle        = sprintf(
    "Species-level · CLR transformation · Top %d species loadings",
    n_arrows
  )
)

print(p_biplot_12)

# ─────────────────────────────────────────────────────────────────────────────
# 5.8.4  Biplot — PC1 vs PC3
# ─────────────────────────────────────────────────────────────────────────────

load_scaled_13 <- get_scaled_loadings(
  rotation  = pca_result$rotation,
  scores_df = pca_df,
  pc_x      = "PC1",
  pc_y      = "PC3",
  n_top     = n_arrows
)

p_biplot_13 <- make_biplot(
  scores_df       = pca_df,
  loadings_scaled = load_scaled_13,
  pc_x_col        = "PC1",
  pc_y_col        = "PC3",
  x_label         = sprintf("PC1 (%.1f%%)", var_exp[1]),
  y_label         = sprintf("PC3 (%.1f%%)", var_exp[3]),
  title           = "PCA Biplot — PC1 vs PC3",
  subtitle        = sprintf(
    "Species-level · CLR transformation · Top %d species loadings",
    n_arrows
  )
)

print(p_biplot_13)

# ─────────────────────────────────────────────────────────────────────────────
# 5.8.5  Combined biplot panel
# ─────────────────────────────────────────────────────────────────────────────

p_biplot_panel <- (p_biplot_12 | p_biplot_13) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "PCA Biplots — Sample Scores & Species Loadings",
    subtitle = sprintf(
      "Arrows = top %d species by PC loading magnitude  ·  CLR-transformed counts",
      n_arrows
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(color = "grey40", size = 10)
    )
  ) &
  theme(legend.position = "bottom")

print(p_biplot_panel)

# ─────────────────────────────────────────────────────────────────────────────
# 5.8.6  Top loadings table — printed for reference
# ─────────────────────────────────────────────────────────────────────────────

loadings_table <- as.data.frame(pca_result$rotation[, 1:3]) %>%
  rownames_to_column("species") %>%
  mutate(
    species_short = sub("^[^|]*\\|", "", species),
    mag_PC1_PC2   = round(sqrt(PC1^2 + PC2^2), 4),
    PC1           = round(PC1, 4),
    PC2           = round(PC2, 4),
    PC3           = round(PC3, 4)
  ) %>%
  arrange(desc(mag_PC1_PC2)) %>%
  select(species_short, PC1, PC2, PC3, mag_PC1_PC2)

cat(sprintf("\n══ Top %d Species Loadings (PC1–PC3) ═══════════════════════════════\n",
            n_arrows))
print(head(loadings_table, n_arrows), row.names = FALSE)
# ─────────────────────────────────────────────────────────────────────────────
# PART 6 — PCA (KO-BASED): CROSS-VALIDATION OF SAMPLE SWAP CHECK
# ─────────────────────────────────────────────────────────────────────────────

# 6.1  Prepare KO count matrix ─────────────────────────────────────────────────

ko_matrix <- ko_counts_final_qc %>%
  filter(!is.na(ko), ko != "") %>%           # remove rows with missing KO IDs
  column_to_rownames("ko") %>%               # set KO IDs as row names
  filter(rowSums(.) > 0)                     # remove all-zero rows

# Transpose → samples as rows, KOs as columns
ko_matrix_t <- t(ko_matrix)

# 6.2  CLR transformation ──────────────────────────────────────────────────────

ko_clr <- ko_matrix_t + 0.5                                       # pseudocount
ko_clr <- sweep(ko_clr, 1, rowSums(ko_clr), "/")                 # relative abundances
ko_clr <- log(ko_clr) - rowMeans(log(ko_clr))                    # CLR

# 6.3  Run PCA ─────────────────────────────────────────────────────────────────

pca_ko_result <- prcomp(ko_clr, center = TRUE, scale. = FALSE)
var_exp_ko    <- summary(pca_ko_result)$importance[2, ] * 100    # % variance explained

# 6.4  Build plotting data frame ───────────────────────────────────────────────

pca_ko_df <- as.data.frame(pca_ko_result$x[, 1:3]) %>%
  rownames_to_column("Sample_ID_FM_Pipeline") %>%
  left_join(
    metadata_final_qc %>% select(Sample_ID_FM_Pipeline, Group),
    by = "Sample_ID_FM_Pipeline"
  )

# 6.5  PC1 vs PC2 ──────────────────────────────────────────────────────────────

p_pca_ko_12 <- ggplot(pca_ko_df,
                      aes(x = PC1, y = PC2,
                          color = Group,
                          label = Sample_ID_FM_Pipeline)) +
  geom_point(size = 3.5, alpha = 0.85) +
  geom_text_repel(
    size          = 3,
    max.overlaps  = 30,
    show.legend   = FALSE,
    segment.color = "grey60"
  ) +
  stat_ellipse(
    aes(group = Group),
    type      = "norm",
    linetype  = "dashed",
    alpha     = 0.5,
    linewidth = 0.7
  ) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title    = "PCA — KO-Based Sample Swap Check (PC1 vs PC2)",
    subtitle = "KEGG Orthology counts · CLR transformation",
    x        = sprintf("PC1 (%.1f%%)", var_exp_ko[1]),
    y        = sprintf("PC2 (%.1f%%)", var_exp_ko[2]),
    color    = "Group"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40"),
    legend.position = "right"
  )

print(p_pca_ko_12)

# 6.6  PC1 vs PC3 ──────────────────────────────────────────────────────────────

p_pca_ko_13 <- ggplot(pca_ko_df,
                      aes(x = PC1, y = PC3,
                          color = Group,
                          label = Sample_ID_FM_Pipeline)) +
  geom_point(size = 3.5, alpha = 0.85) +
  geom_text_repel(
    size          = 3,
    max.overlaps  = 30,
    show.legend   = FALSE,
    segment.color = "grey60"
  ) +
  stat_ellipse(
    aes(group = Group),
    type      = "norm",
    linetype  = "dashed",
    alpha     = 0.5,
    linewidth = 0.7
  ) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title    = "PCA — KO-Based Sample Swap Check (PC1 vs PC3)",
    subtitle = "KEGG Orthology counts · CLR transformation",
    x        = sprintf("PC1 (%.1f%%)", var_exp_ko[1]),
    y        = sprintf("PC3 (%.1f%%)", var_exp_ko[3]),
    color    = "Group"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40"),
    legend.position = "right"
  )

print(p_pca_ko_13)

# 6.7  Scree plot ──────────────────────────────────────────────────────────────

scree_ko_df <- data.frame(
  PC      = seq_along(var_exp_ko),
  Var_Pct = var_exp_ko
)

p_scree_ko <- ggplot(scree_ko_df[1:15, ], aes(x = PC, y = Var_Pct)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_line(color = "darkred", linewidth = 0.8) +
  geom_point(color = "darkred", size = 2) +
  labs(
    title = "Scree Plot — KO PCA Variance Explained per PC",
    x     = "Principal Component",
    y     = "Variance Explained (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

print(p_scree_ko)

# ─────────────────────────────────────────────────────────────────────────────
# 6.8  SIDE-BY-SIDE COMPARISON — Taxonomic vs KO (PC1 vs PC2)
# ─────────────────────────────────────────────────────────────────────────────

library(patchwork)

# Add a source label to each PCA data frame
pca_combined <- bind_rows(
  pca_df    %>% mutate(Source = "Taxonomic (Species)"),
  pca_ko_df %>% mutate(Source = "Functional (KO)")
)

# Normalise PC axes direction isn't guaranteed between runs,
# so plot them separately but display together with patchwork

p_compare <- (p_pca_12 + labs(title = "Taxonomic PCA")) +
  (p_pca_ko_12 + labs(title = "Functional (KO) PCA")) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(p_compare)


# ------------------------------------------------- Alpha --------------
# ─────────────────────────────────────────────────────────────────────────────
# PART 7 — ALPHA DIVERSITY: RAREFACTION + METRICS + PLOTS + STATS
# ─────────────────────────────────────────────────────────────────────────────

# 7.1  Build phyloseq object (species-level) ───────────────────────────────────

sample_cols_tax <- setdiff(colnames(tax_counts_final_qc), tax_cols)

species_df_alpha <- tax_counts_final_qc %>%
  filter(!is.na(s), s != "") %>%
  mutate(tax_id = make.unique(paste(g, s, sep = "|")))

# OTU table — taxa × samples
otu_mat_alpha <- species_df_alpha %>%
  select(tax_id, all_of(sample_cols_tax)) %>%
  column_to_rownames("tax_id") %>%
  as.matrix()

# Taxonomy table — taxa × ranks
tax_mat_alpha <- species_df_alpha %>%
  select(tax_id, all_of(tax_cols)) %>%
  column_to_rownames("tax_id") %>%
  as.matrix()

# Sample metadata
samp_df_alpha <- metadata_final_qc %>%
  column_to_rownames("Sample_ID_FM_Pipeline") %>%
  as.data.frame()

# Assemble phyloseq
ps_alpha <- phyloseq(
  otu_table(otu_mat_alpha, taxa_are_rows = TRUE),
  tax_table(tax_mat_alpha),
  sample_data(samp_df_alpha)
)

cat(sprintf("\nPhyloseq built: %d taxa × %d samples\n",
            ntaxa(ps_alpha), nsamples(ps_alpha)))

# ─────────────────────────────────────────────────────────────────────────────
# 7.2  Sequencing depth check — before rarefaction
# ─────────────────────────────────────────────────────────────────────────────
# NOTE: phyloseq applies make.names() to numeric sample IDs
#       (e.g. "12345" → "X12345") — fix_id() reverses this before joining

depth_df <- data.frame(
  Sample_ID_FM_Pipeline = fix_id(sample_names(ps_alpha)),
  Depth                 = sample_sums(ps_alpha)
) %>%
  left_join(
    metadata_final_qc %>% select(Sample_ID_FM_Pipeline, Group),
    by = "Sample_ID_FM_Pipeline"
  ) %>%
  arrange(Depth)

# Validate join
if (any(is.na(depth_df$Group))) {
  warning("Depth join failed for some samples — check sample IDs!")
  print(depth_df %>% filter(is.na(Group)))
} else {
  cat("✔ Depth join OK — all samples matched to metadata\n")
}

rare_depth <- min(depth_df$Depth)

cat(sprintf(
  "Rarefaction depth (minimum): %d  |  Max: %d  |  Median: %.0f\n",
  rare_depth, max(depth_df$Depth), median(depth_df$Depth)
))

p_depth <- ggplot(depth_df,
                  aes(x = reorder(Sample_ID_FM_Pipeline, Depth),
                      y = Depth,
                      fill = Group)) +
  geom_col(alpha = 0.85) +
  geom_hline(yintercept = rare_depth,
             color = "darkred", linetype = "dashed", linewidth = 0.8) +
  annotate("text",
           x     = 1,
           y     = rare_depth,
           label = sprintf("Rarefaction depth: %d", rare_depth),
           vjust = -0.5, hjust = 0, color = "darkred", size = 3.5) +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title    = "Sequencing Depth per Sample",
    subtitle = "Red dashed line = rarefaction threshold (minimum depth)",
    x        = "Sample",
    y        = "Read Count",
    fill     = "Group"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40"),
    axis.text.x   = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

print(p_depth)

# ─────────────────────────────────────────────────────────────────────────────
# 7.3  Rarefaction — phyloseq::rarefy_even_depth
# ─────────────────────────────────────────────────────────────────────────────
# NOTE: rngseed = 42 used directly — cleaner than set.seed() + rngseed = FALSE

ps_rare <- rarefy_even_depth(
  ps_alpha,
  sample.size = rare_depth,
  rngseed     = 42,
  replace     = FALSE,
  trimOTUs    = TRUE,
  verbose     = TRUE
)

cat(sprintf(
  "\nAfter rarefaction: %d taxa × %d samples | Depth: %d reads/sample\n",
  ntaxa(ps_rare), nsamples(ps_rare), rare_depth
))

# ─────────────────────────────────────────────────────────────────────────────
# 7.4  Estimate alpha diversity metrics
# ─────────────────────────────────────────────────────────────────────────────
# NOTE: fix_id() applied to rownames — phyloseq's make.names() adds "X" prefix
#       to numeric IDs, which breaks the left_join with metadata if not reversed

alpha_div <- estimate_richness(
  ps_rare,
  measures = c("Observed", "Chao1", "Shannon", "Simpson")
) %>%
  rownames_to_column("Sample_ID_FM_Pipeline") %>%
  mutate(Sample_ID_FM_Pipeline = fix_id(Sample_ID_FM_Pipeline)) %>%
  left_join(
    metadata_final_qc %>% select(Sample_ID_FM_Pipeline, Group),
    by = "Sample_ID_FM_Pipeline"
  )

# Validate join
if (any(is.na(alpha_div$Group))) {
  warning("Alpha div join failed — some samples missing Group assignment!")
  print(alpha_div %>% filter(is.na(Group)) %>% select(Sample_ID_FM_Pipeline, Group))
} else {
  cat(sprintf("✔ Alpha diversity computed for %d samples across %d groups\n",
              nrow(alpha_div), n_distinct(alpha_div$Group)))
}

# Long format for plotting
alpha_long <- alpha_div %>%
  select(Sample_ID_FM_Pipeline, Group, Observed, Chao1, Shannon, Simpson) %>%
  pivot_longer(
    cols      = c(Observed, Chao1, Shannon, Simpson),
    names_to  = "Metric",
    values_to = "Value"
  ) %>%
  mutate(Metric = factor(Metric,
                         levels = c("Observed", "Chao1", "Shannon", "Simpson")))

# ─────────────────────────────────────────────────────────────────────────────
# 7.5  Colour palette & comparison setup
# ─────────────────────────────────────────────────────────────────────────────
# NOTE: 2 groups only — list() used directly instead of combn()

groups_ordered   <- sort(unique(alpha_long$Group))
group_colors     <- setNames(
  RColorBrewer::brewer.pal(3, "Set1")[1:2],
  groups_ordered
)
comparisons_list <- list(groups_ordered)   # single pair for 2 groups

cat(sprintf("\nGroups: %s\n", paste(groups_ordered, collapse = " vs ")))

# ─────────────────────────────────────────────────────────────────────────────
# 7.6  Faceted violin + boxplot — all 4 metrics
# ─────────────────────────────────────────────────────────────────────────────
# NOTE: stat_compare_means(kruskal.test) REMOVED — redundant for 2 groups
#       (K-W equals Wilcoxon) and causes label misplacement in free_y facets

p_alpha_all <- ggplot(alpha_long,
                      aes(x     = Group,
                          y     = Value,
                          fill  = Group,
                          color = Group)) +
  geom_violin(alpha = 0.35, linewidth = 0.6, trim = FALSE) +
  geom_boxplot(width         = 0.18,
               alpha         = 0.80,
               outlier.shape = NA,
               color         = "black") +
  geom_jitter(width = 0.08, size = 1.8, alpha = 0.70, show.legend = FALSE) +
  stat_compare_means(
    method       = "wilcox.test",
    comparisons  = comparisons_list,
    label        = "p.format",
    tip.length   = 0.01,
    bracket.size = 0.4,
    size         = 3.8
  ) +
  scale_fill_manual(values  = group_colors) +
  scale_color_manual(values = group_colors) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 2) +
  labs(
    title    = "Alpha Diversity — Rarefied Data",
    subtitle = sprintf(
      "Rarefied to %d reads/sample  ·  Wilcoxon rank-sum p-value shown",
      rare_depth),
    x     = NULL,
    y     = "Diversity Value",
    fill  = "Group",
    color = "Group"
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

print(p_alpha_all)

# ─────────────────────────────────────────────────────────────────────────────
# 7.7  Individual metric plots — 2×2 publication panel
# ─────────────────────────────────────────────────────────────────────────────

plot_alpha_metric <- function(metric_name, y_label) {
  
  dat <- alpha_long %>% filter(Metric == metric_name)
  
  ggplot(dat, aes(x = Group, y = Value, fill = Group, color = Group)) +
    geom_violin(alpha = 0.30, linewidth = 0.6, trim = FALSE) +
    geom_boxplot(width         = 0.20,
                 alpha         = 0.80,
                 outlier.shape = NA,
                 color         = "black") +
    geom_jitter(width = 0.08, size = 2.2, alpha = 0.75, show.legend = FALSE) +
    stat_compare_means(
      method      = "wilcox.test",
      comparisons = comparisons_list,
      label       = "p.format",
      tip.length  = 0.01,
      size        = 3.8
    ) +
    scale_fill_manual(values  = group_colors) +
    scale_color_manual(values = group_colors) +
    labs(title = metric_name, x = NULL, y = y_label) +
    theme_bw(base_size = 13) +
    theme(
      plot.title      = element_text(face = "bold", hjust = 0.5),
      legend.position = "none",
      axis.text.x     = element_text(angle = 20, hjust = 1)
    )
}

p_observed <- plot_alpha_metric("Observed", "Observed Species (Richness)")
p_chao1    <- plot_alpha_metric("Chao1",    "Chao1 Estimated Richness")
p_shannon  <- plot_alpha_metric("Shannon",  "Shannon Diversity Index")
p_simpson  <- plot_alpha_metric("Simpson",  "Simpson Diversity Index")

p_alpha_panel <- (p_observed | p_chao1) / (p_shannon | p_simpson) +
  plot_annotation(
    title    = "Alpha Diversity — Individual Metrics",
    subtitle = sprintf("Rarefied to %d reads/sample · Wilcoxon p-value shown",
                       rare_depth),
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(color = "grey40", size = 10)
    )
  )

print(p_alpha_panel)

# ─────────────────────────────────────────────────────────────────────────────
# 7.8  Statistical tests — full results tables
# ─────────────────────────────────────────────────────────────────────────────

metrics_list <- c("Observed", "Chao1", "Shannon", "Simpson")
grp_a        <- groups_ordered[1]
grp_b        <- groups_ordered[2]

# ── Kruskal-Wallis ─────────────────────────────────────────────────────────────

kw_results <- map_dfr(metrics_list, function(m) {
  dat  <- alpha_div %>% filter(!is.na(.data[[m]]))
  test <- kruskal.test(reformulate("Group", response = m), data = dat)
  tibble(
    Metric    = m,
    Test      = "Kruskal-Wallis",
    Statistic = round(test$statistic, 4),
    df        = test$parameter,
    p_value   = round(test$p.value, 4),
    Sig       = case_when(
      test$p.value < 0.001 ~ "***",
      test$p.value < 0.01  ~ "**",
      test$p.value < 0.05  ~ "*",
      TRUE                 ~ "ns"
    )
  )
})

# ── Wilcoxon rank-sum ──────────────────────────────────────────────────────────
# NOTE: wilcox.test() used directly — cleaner for 2 groups, returns W statistic
# NOTE: exact = FALSE required — microbiome data always has ties

wx_results <- map_dfr(metrics_list, function(m) {
  dat    <- alpha_div %>% filter(!is.na(.data[[m]]))
  vals_a <- dat %>% filter(Group == grp_a) %>% pull(.data[[m]])
  vals_b <- dat %>% filter(Group == grp_b) %>% pull(.data[[m]])
  test   <- wilcox.test(vals_a, vals_b, exact = FALSE)
  tibble(
    Metric     = m,
    Test       = "Wilcoxon rank-sum",
    Comparison = paste(grp_a, "vs", grp_b),
    W          = round(test$statistic, 4),
    p_value    = round(test$p.value, 4),
    Sig        = case_when(
      test$p.value < 0.001 ~ "***",
      test$p.value < 0.01  ~ "**",
      test$p.value < 0.05  ~ "*",
      TRUE                 ~ "ns"
    )
  )
})

cat("\n══ Kruskal-Wallis Results (note: equivalent to Wilcoxon for 2 groups) ══\n")
print(kw_results, n = Inf)

cat("\n══ Wilcoxon Rank-Sum Results ════════════════════════════════════════════\n")
print(wx_results, n = Inf)

# ── Descriptive statistics per group ──────────────────────────────────────────

alpha_summary <- alpha_long %>%
  group_by(Group, Metric) %>%
  summarise(
    N      = n(),
    Mean   = round(mean(Value,   na.rm = TRUE), 4),
    Median = round(median(Value, na.rm = TRUE), 4),
    SD     = round(sd(Value,     na.rm = TRUE), 4),
    Min    = round(min(Value,    na.rm = TRUE), 4),
    Max    = round(max(Value,    na.rm = TRUE), 4),
    .groups = "drop"
  )

cat("\n══ Descriptive Statistics per Group ════════════════════════════════════\n")
print(alpha_summary, n = Inf)


#------------------------------------------ BETA -----------------------------

# ─────────────────────────────────────────────────────────────────────────────
# PART 8 — BETA DIVERSITY: DISTANCES + ORDINATION + HEATMAP + STATISTICS
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# 8.1  Prepare rarefied OTU matrix from ps_rare (Part 7)
# ─────────────────────────────────────────────────────────────────────────────
# NOTE: fix_id() reverses phyloseq make.names() on numeric sample IDs

otu_rare_mat <- as.data.frame(otu_table(ps_rare))
if (taxa_are_rows(ps_rare)) otu_rare_mat <- t(otu_rare_mat)
otu_rare_mat <- as.matrix(otu_rare_mat)
rownames(otu_rare_mat) <- fix_id(rownames(otu_rare_mat))

# Relative abundance matrix
otu_rel_mat <- sweep(otu_rare_mat, 1, rowSums(otu_rare_mat), "/")

# Align metadata to OTU matrix row order
meta_beta <- metadata_final_qc %>%
  filter(Sample_ID_FM_Pipeline %in% rownames(otu_rare_mat)) %>%
  arrange(match(Sample_ID_FM_Pipeline, rownames(otu_rare_mat)))

stopifnot(all(rownames(otu_rare_mat) == meta_beta$Sample_ID_FM_Pipeline))
cat(sprintf("✔ OTU matrix aligned: %d samples × %d species\n",
            nrow(otu_rare_mat), ncol(otu_rare_mat)))

# ─────────────────────────────────────────────────────────────────────────────
# 8.2  Compute distance matrices
# ─────────────────────────────────────────────────────────────────────────────

# Bray-Curtis — abundance-based
dist_bray <- vegdist(otu_rel_mat, method = "bray")

# Jaccard — presence/absence
dist_jacc <- vegdist(otu_rare_mat, method = "jaccard", binary = TRUE)

# Weighted UniFrac — requires phylogenetic tree in ps_rare
# For shotgun metagenomics: add a tree with phy_tree(ps_rare) <- your_tree
# For 16S data: tree is typically already present in the phyloseq object
unifrac_available <- !is.null(tryCatch(phy_tree(ps_rare), error = function(e) NULL))

if (unifrac_available) {
  dist_unifrac_w <- UniFrac(ps_rare, weighted = TRUE, normalized = TRUE)
  cat("✔ Weighted UniFrac computed\n")
} else {
  cat("⚠ Weighted UniFrac skipped — no phylogenetic tree found in ps_rare\n")
  cat("  To enable: add a tree object → phy_tree(ps_rare) <- read.tree('your_tree.nwk')\n")
}

# Build distance list — UniFrac added only if tree is available
dist_list <- list(
  "Bray-Curtis"      = dist_bray,
  "Jaccard"          = dist_jacc
)
if (unifrac_available) {
  dist_list[["Weighted UniFrac"]] <- dist_unifrac_w
}

cat(sprintf("\nDistance matrices ready: %s\n",
            paste(names(dist_list), collapse = ", ")))

# ─────────────────────────────────────────────────────────────────────────────
# 8.3  Pre-compute PERMANOVA + ANOSIM for ALL distances
#      *** Must run BEFORE PCoA plots so stats can be annotated on each plot ***
# ─────────────────────────────────────────────────────────────────────────────

sig_stars <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

stats_by_dist <- map(names(dist_list), function(dist_name) {
  
  # PERMANOVA
  set.seed(42)
  perm <- adonis2(
    dist_list[[dist_name]] ~ Group,
    data         = meta_beta,
    permutations = 999,
    by           = "margin"
  )
  perm_R2 <- round(perm$R2[1],          3)
  perm_p  <- perm$`Pr(>F)`[1]
  perm_F  <- round(perm$F[1],           3)
  
  # ANOSIM
  set.seed(42)
  ano    <- anosim(dist_list[[dist_name]],
                   grouping     = meta_beta$Group,
                   permutations = 999)
  ano_R  <- round(ano$statistic, 3)
  ano_p  <- ano$signif
  
  cat(sprintf(
    "[%s] PERMANOVA: R²=%.3f, F=%.3f, p=%.4f%s | ANOSIM: R=%.3f, p=%.4f%s\n",
    dist_name,
    perm_R2, perm_F, perm_p, sig_stars(perm_p),
    ano_R, ano_p, sig_stars(ano_p)
  ))
  
  list(
    permanova_R2  = perm_R2,
    permanova_F   = perm_F,
    permanova_p   = perm_p,
    permanova_sig = sig_stars(perm_p),
    anosim_R      = ano_R,
    anosim_p      = ano_p,
    anosim_sig    = sig_stars(ano_p),
    # Formatted annotation string for plot
    annot_text    = sprintf(
      "PERMANOVA: R\u00B2 = %.3f, p = %.4f %s\nANOSIM:    R = %.3f,   p = %.4f %s",
      perm_R2, perm_p, sig_stars(perm_p),
      ano_R,   ano_p,  sig_stars(ano_p)
    )
  )
}) %>% setNames(names(dist_list))

# ─────────────────────────────────────────────────────────────────────────────
# 8.4  Helper — ordination plot with optional stats annotation box
# ─────────────────────────────────────────────────────────────────────────────

plot_ord_stats <- function(scores_df, x_col, y_col,
                           x_label, y_label,
                           title, subtitle,
                           stats_text = NULL) {
  
  p <- ggplot(scores_df,
              aes(x     = .data[[x_col]],
                  y     = .data[[y_col]],
                  color = Group,
                  label = Sample_ID_FM_Pipeline)) +
    geom_point(size = 3.5, alpha = 0.85) +
    geom_text_repel(
      size          = 2.8,
      max.overlaps  = 30,
      show.legend   = FALSE,
      segment.color = "grey60"
    ) +
    stat_ellipse(
      aes(group = Group),
      type      = "norm",
      linetype  = "dashed",
      alpha     = 0.6,
      linewidth = 0.7
    ) +
    scale_color_manual(values = group_colors) +
    labs(
      title    = title,
      subtitle = subtitle,
      x        = x_label,
      y        = y_label,
      color    = "Group"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title      = element_text(face = "bold"),
      plot.subtitle   = element_text(color = "grey40", size = 10),
      legend.position = "right"
    )
  
  # Add PERMANOVA + ANOSIM annotation box (bottom-right corner)
  if (!is.null(stats_text)) {
    p <- p +
      annotate(
        "label",
        x          = Inf,
        y          = -Inf,
        label      = stats_text,
        hjust      = 1.03,
        vjust      = -0.15,
        size       = 3.0,
        color      = "black",
        fill       = "white",
        alpha      = 0.88,
        label.size = 0.35,
        family     = "mono",
        fontface   = "plain"
      )
  }
  
  p
}

# ─────────────────────────────────────────────────────────────────────────────
# 8.5  PCoA — with PERMANOVA + ANOSIM annotations
# ─────────────────────────────────────────────────────────────────────────────

pcoa_plots <- list()

for (dist_name in names(dist_list)) {
  
  # Run PCoA
  pcoa_res <- pcoa(dist_list[[dist_name]])
  
  # Variance explained (positive eigenvalues only)
  eig      <- pcoa_res$values$Eigenvalues
  var_pcoa <- round(eig[eig > 0] / sum(eig[eig > 0]) * 100, 1)
  
  # Build scores data frame
  pcoa_df <- as.data.frame(pcoa_res$vectors[, 1:2]) %>%
    setNames(c("Axis1", "Axis2")) %>%
    rownames_to_column("Sample_ID_FM_Pipeline") %>%
    mutate(Sample_ID_FM_Pipeline = fix_id(Sample_ID_FM_Pipeline)) %>%
    left_join(
      metadata_final_qc %>% select(Sample_ID_FM_Pipeline, Group),
      by = "Sample_ID_FM_Pipeline"
    )
  
  # Retrieve pre-computed stats annotation
  annot <- stats_by_dist[[dist_name]]$annot_text
  
  p <- plot_ord_stats(
    scores_df  = pcoa_df,
    x_col      = "Axis1",
    y_col      = "Axis2",
    x_label    = sprintf("PCoA Axis 1 (%.1f%%)", var_pcoa[1]),
    y_label    = sprintf("PCoA Axis 2 (%.1f%%)", var_pcoa[2]),
    title      = sprintf("PCoA — %s Distance", dist_name),
    subtitle   = "Rarefied species counts | 999 permutations",
    stats_text = annot
  )
  
  pcoa_plots[[dist_name]] <- p
  print(p)
}

# Combined panel — 2 columns (3 plots if UniFrac available)
n_cols_pcoa <- ifelse(length(pcoa_plots) == 3, 3, 2)

p_pcoa_panel <- wrap_plots(pcoa_plots, ncol = n_cols_pcoa) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "PCoA — Beta Diversity (Rarefied Data)",
    subtitle = "",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(color = "grey40", size = 10)
    )
  ) &
  theme(legend.position = "bottom")

print(p_pcoa_panel)

# ─────────────────────────────────────────────────────────────────────────────
# 8.6  NMDS — Non-metric Multidimensional Scaling
# ─────────────────────────────────────────────────────────────────────────────

nmds_plots  <- list()
nmds_stress <- list()

for (dist_name in names(dist_list)) {
  
  set.seed(42)
  nmds_res <- metaMDS(
    dist_list[[dist_name]],
    k             = 2,
    trymax        = 100,
    trace         = FALSE,
    autotransform = FALSE
  )
  
  stress_val   <- nmds_res$stress
  stress_label <- case_when(
    stress_val < 0.10 ~ "excellent",
    stress_val < 0.20 ~ "good",
    stress_val < 0.30 ~ "acceptable",
    TRUE              ~ "poor — interpret with caution"
  )
  
  nmds_stress[[dist_name]] <- stress_val
  cat(sprintf("\nNMDS [%s]: stress = %.4f (%s)\n",
              dist_name, stress_val, stress_label))
  
  nmds_df <- as.data.frame(scores(nmds_res, display = "sites")) %>%
    rownames_to_column("Sample_ID_FM_Pipeline") %>%
    mutate(Sample_ID_FM_Pipeline = fix_id(Sample_ID_FM_Pipeline)) %>%
    left_join(
      metadata_final_qc %>% select(Sample_ID_FM_Pipeline, Group),
      by = "Sample_ID_FM_Pipeline"
    )
  
  p <- plot_ord_stats(
    scores_df  = nmds_df,
    x_col      = "NMDS1",
    y_col      = "NMDS2",
    x_label    = "NMDS1",
    y_label    = "NMDS2",
    title      = sprintf("NMDS — %s Distance", dist_name),
    subtitle   = sprintf("Stress = %.4f (%s) | Rarefied species counts",
                         stress_val, stress_label),
    stats_text = NULL   # stats shown on PCoA; NMDS keeps stress annotation only
  ) +
    annotate("text",
             x     = -Inf,
             y     = Inf,
             label = sprintf("Stress: %.4f", stress_val),
             hjust = -0.1,
             vjust = 1.4,
             size  = 3.5,
             color = "grey30",
             fontface = "italic")
  
  nmds_plots[[dist_name]] <- p
  print(p)
}

# Combined NMDS panel
p_nmds_panel <- wrap_plots(nmds_plots, ncol = n_cols_pcoa) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "NMDS — Beta Diversity (Rarefied Data)",
    theme = theme(plot.title = element_text(face = "bold", size = 15))
  ) &
  theme(legend.position = "bottom")

print(p_nmds_panel)

# ─────────────────────────────────────────────────────────────────────────────
# 8.7  ComplexHeatmap — Top variable species
# ─────────────────────────────────────────────────────────────────────────────

n_top <- 50

top_species <- apply(otu_rel_mat, 2, var) %>%
  sort(decreasing = TRUE) %>%
  head(n_top) %>%
  names()

heat_mat   <- t(otu_rel_mat[, top_species])
heat_mat_z <- t(scale(t(heat_mat)))
heat_mat_z <- pmax(pmin(heat_mat_z, 3), -3)

rownames(heat_mat_z) <- sub("^[^|]*\\|", "", rownames(heat_mat_z))
colnames(heat_mat_z) <- fix_id(colnames(heat_mat_z))

hclust_cols <- hclust(dist_bray,        method = "ward.D2")
hclust_rows <- hclust(dist(heat_mat_z), method = "ward.D2")

col_fun <- colorRamp2(
  c(-3, -1.5, 0, 1.5, 3),
  c("#2166AC", "#92C5DE", "white", "#F4A582", "#D6604D")
)

col_anno <- HeatmapAnnotation(
  Group = meta_beta$Group,
  col   = list(Group = group_colors),
  annotation_name_gp      = gpar(fontsize = 11, fontface = "bold"),
  annotation_legend_param = list(
    Group = list(title_gp = gpar(fontface = "bold"))
  )
)

ht <- Heatmap(
  heat_mat_z,
  name                 = "Z-score",
  col                  = col_fun,
  top_annotation       = col_anno,
  cluster_rows         = hclust_rows,
  cluster_columns      = hclust_cols,
  show_column_names    = TRUE,
  column_names_gp      = gpar(fontsize = 8),
  show_row_names       = TRUE,
  row_names_gp         = gpar(fontsize = 7.5),
  row_names_side       = "left",
  row_dend_side        = "right",
  column_title         = sprintf(
    "Beta Diversity Heatmap — Top %d Variable Species (Z-score, Bray-Curtis clustering)",
    n_top
  ),
  column_title_gp      = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    title    = "Z-score",
    at       = c(-3, -1.5, 0, 1.5, 3),
    title_gp = gpar(fontface = "bold")
  ),
  border = TRUE
)

draw(ht, merge_legend = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# 8.8  Homogeneity of dispersion — betadisper (PERMANOVA assumption check)
# ─────────────────────────────────────────────────────────────────────────────

betadisp_results <- map_dfr(names(dist_list), function(dist_name) {
  
  beta_d <- betadisper(dist_list[[dist_name]], group = meta_beta$Group)
  set.seed(42)
  perm_d <- permutest(beta_d, permutations = 999)
  
  tibble(
    Distance = dist_name,
    F_stat   = round(perm_d$tab[1, "F"],        4),
    p_value  = round(perm_d$tab[1, "Pr(>F)"],   4),
    Result   = case_when(
      perm_d$tab[1, "Pr(>F)"] < 0.05 ~
        "⚠ Unequal dispersion — PERMANOVA result may reflect spread, not location",
      TRUE ~
        "✔ Homogeneous dispersion — PERMANOVA result reliable"
    )
  )
})

# ─────────────────────────────────────────────────────────────────────────────
# 8.9  Build full results tables from pre-computed stats (8.3)
# ─────────────────────────────────────────────────────────────────────────────

permanova_results <- map_dfr(names(stats_by_dist), function(dist_name) {
  s <- stats_by_dist[[dist_name]]
  tibble(
    Distance = dist_name,
    R2       = s$permanova_R2,
    F_stat   = s$permanova_F,
    p_value  = s$permanova_p,
    Sig      = s$permanova_sig
  )
})

anosim_results <- map_dfr(names(stats_by_dist), function(dist_name) {
  s <- stats_by_dist[[dist_name]]
  tibble(
    Distance = dist_name,
    R_stat   = s$anosim_R,
    p_value  = s$anosim_p,
    Strength = case_when(
      s$anosim_R > 0.75 ~ "strong",
      s$anosim_R > 0.50 ~ "moderate",
      s$anosim_R > 0.25 ~ "weak",
      TRUE              ~ "negligible"
    ),
    Sig = s$anosim_sig
  )
})

# ─────────────────────────────────────────────────────────────────────────────
# 8.10  Print full statistical summary
# ─────────────────────────────────────────────────────────────────────────────

cat("\n══════════════════════════════════════════════════════════════════════════\n")
cat("  BETA DIVERSITY — COMPLETE STATISTICAL SUMMARY\n")
cat("══════════════════════════════════════════════════════════════════════════\n")

cat("\n── PERMANOVA (adonis2) ──────────────────────────────────────────────────\n")
cat("Tests: do groups differ in community composition? (centroid differences)\n\n")
print(permanova_results, n = Inf)

cat("\n── ANOSIM ───────────────────────────────────────────────────────────────\n")
cat("Tests: between-group vs within-group dissimilarity (R: -1 to 1)\n\n")
print(anosim_results, n = Inf)

cat("\n── Homogeneity of Dispersion (betadisper) ───────────────────────────────\n")
cat("PERMANOVA assumption: groups must have equal within-group spread\n\n")
print(betadisp_results, n = Inf)


# ─────────────────────────────────────────────────────────────────────────────
# PART 9 — DIFFERENTIAL ABUNDANCE: MaAsLin2 + VOLCANO PLOT
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# 9.1  Prepare MaAsLin2 input
#      Uses unrarefied filtered counts — MaAsLin2 handles normalisation internally
#      Rarefied data (ps_rare) is NOT recommended for differential abundance
# ─────────────────────────────────────────────────────────────────────────────

sample_cols_tax <- setdiff(colnames(tax_counts_final_qc), tax_cols)

species_maaslin <- tax_counts_final_qc %>%
  filter(!is.na(s), s != "") %>%
  mutate(tax_id = make.unique(paste(g, s, sep = "|")))

# Features matrix: samples × species  (MaAsLin2 expects rows = samples)
maaslin_features <- species_maaslin %>%
  select(tax_id, all_of(sample_cols_tax)) %>%
  column_to_rownames("tax_id") %>%
  t() %>%
  as.data.frame()

rownames(maaslin_features) <- fix_id(rownames(maaslin_features))

# Metadata: samples × variables
maaslin_meta <- metadata_final_qc %>%
  column_to_rownames("Sample_ID_FM_Pipeline") %>%
  as.data.frame()

# Align to common samples
common_s       <- intersect(rownames(maaslin_features), rownames(maaslin_meta))
maaslin_features <- maaslin_features[common_s, ]
maaslin_meta     <- maaslin_meta[common_s, , drop = FALSE]

cat(sprintf("\nMaAsLin2 input: %d samples × %d features\n",
            nrow(maaslin_features), ncol(maaslin_features)))

# Set reference group
# Positive coefficient  → higher in comparison group (Symphiome)
# Negative coefficient  → higher in reference group  (Cp Challenge Control)
ref_group        <- "Cp Challenge Control"
comparison_group <- setdiff(unique(maaslin_meta$Group), ref_group)
maaslin_meta$Group <- relevel(factor(maaslin_meta$Group), ref = ref_group)

cat(sprintf("Reference : %s\nComparison: %s\n", ref_group, comparison_group))

# ─────────────────────────────────────────────────────────────────────────────
# 9.2  Run MaAsLin2
# ─────────────────────────────────────────────────────────────────────────────

maaslin_output_dir <- "maaslin2_output"

set.seed(42)
maaslin_fit <- Maaslin2(
  input_data      = maaslin_features,
  input_metadata  = maaslin_meta,
  output          = maaslin_output_dir,
  fixed_effects   = "Group",
  normalization   = "TSS",     # Total Sum Scaling → relative abundance
  transform       = "LOG",     # Log transformation
  analysis_method = "LM",      # Linear model
  min_abundance   = 0,         # Already filtered in Part 4.7
  min_prevalence  = 0,         # Already filtered in Part 4.7
  max_significance = 0.25,     # MaAsLin2 default q-value threshold
  correction      = "BH",      # Benjamini-Hochberg FDR
  standardize     = FALSE,
  plot_heatmap    = FALSE,      # Custom plots below
  plot_scatter    = FALSE,
  cores           = 1
)

cat(sprintf("\n✔ MaAsLin2 complete. Results saved to: %s/\n", maaslin_output_dir))

# ─────────────────────────────────────────────────────────────────────────────
# 9.3  Process results
# ─────────────────────────────────────────────────────────────────────────────

# Significance thresholds (adjustable)
q_threshold    <- 0.25    # MaAsLin2 standard; use 0.05 for stricter analysis
coef_threshold <- 0.5     # minimum effect size for volcano dashed lines

maaslin_res <- maaslin_fit$results %>%
  as_tibble() %>%
  filter(metadata == "Group") %>%
  arrange(qval, pval) %>%
  mutate(
    # Clean display name: keep only the species part after "|"
    feature_clean = sub("^[^|]*\\|", "", feature),
    # Guard against qval = 0 or NA before -log10
    qval_safe     = pmax(qval, 1e-10),
    neg_log10_q   = -log10(qval_safe),
    # Direction labels
    Direction = case_when(
      qval < q_threshold & coef > 0 ~ paste("Enriched in", comparison_group),
      qval < q_threshold & coef < 0 ~ paste("Enriched in", ref_group),
      TRUE                           ~ "Not significant"
    )
  )

# Summary counts
n_sig_25  <- sum(maaslin_res$qval < 0.25, na.rm = TRUE)
n_sig_05  <- sum(maaslin_res$qval < 0.05, na.rm = TRUE)
n_enr_comp <- sum(maaslin_res$qval < q_threshold & maaslin_res$coef > 0, na.rm = TRUE)
n_enr_ref  <- sum(maaslin_res$qval < q_threshold & maaslin_res$coef < 0, na.rm = TRUE)

cat(sprintf(
  "\nFeatures tested : %d\nSignificant q<0.25 : %d\nSignificant q<0.05 : %d\n",
  nrow(maaslin_res), n_sig_25, n_sig_05
))
cat(sprintf(
  "Enriched in %-30s: %d\nEnriched in %-30s: %d\n",
  comparison_group, n_enr_comp,
  ref_group,        n_enr_ref
))

# Print top hits
cat("\n══ Top Significant Features ═════════════════════════════════════════════\n")
print(
  maaslin_res %>%
    filter(qval < q_threshold) %>%
    select(feature_clean, coef, stderr, pval, qval, Direction) %>%
    arrange(qval),
  n = 30
)

# ─────────────────────────────────────────────────────────────────────────────
# 9.4  Colour scheme (inherits group_colors from Part 7)
# ─────────────────────────────────────────────────────────────────────────────

volcano_colors <- c(
  setNames(group_colors[comparison_group], paste("Enriched in", comparison_group)),
  setNames(group_colors[ref_group],        paste("Enriched in", ref_group)),
  "Not significant" = "grey72"
)

# Top labels: up to 12 per direction, prioritised by lowest q then largest |coef|
top_labels <- maaslin_res %>%
  filter(qval < q_threshold) %>%
  group_by(Direction) %>%
  arrange(qval, desc(abs(coef))) %>%
  slice_head(n = 12) %>%
  ungroup()

# ─────────────────────────────────────────────────────────────────────────────
# 9.5  Volcano plot
# ─────────────────────────────────────────────────────────────────────────────

p_volcano <- ggplot(maaslin_res,
                    aes(x     = coef,
                        y     = neg_log10_q,
                        color = Direction,
                        size  = abs(coef))) +
  # Background points
  geom_point(alpha = 0.72) +
  # Labels for significant features
  geom_text_repel(
    data          = top_labels,
    aes(label     = feature_clean),
    size          = 2.9,
    max.overlaps  = 30,
    segment.color = "grey50",
    segment.size  = 0.3,
    box.padding   = 0.45,
    show.legend   = FALSE
  ) +
  # Horizontal q-value threshold line
  geom_hline(
    yintercept = -log10(q_threshold),
    linetype   = "dashed",
    color      = "grey35",
    linewidth  = 0.6
  ) +
  # Vertical effect size threshold lines
  geom_vline(
    xintercept = c(-coef_threshold, coef_threshold),
    linetype   = "dashed",
    color      = "grey35",
    linewidth  = 0.6
  ) +
  # Threshold annotations
  annotate("text",
           x = max(maaslin_res$coef, na.rm = TRUE),
           y = -log10(q_threshold),
           label    = sprintf("q = %.2f", q_threshold),
           hjust    = 1, vjust = -0.4,
           size     = 3.2, color = "grey35", fontface = "italic") +
  # Enrichment count labels in corners
  annotate("text",
           x = min(maaslin_res$coef, na.rm = TRUE),
           y = max(maaslin_res$neg_log10_q, na.rm = TRUE),
           label    = sprintf("n = %d", n_enr_ref),
           hjust    = 0, vjust = 1,
           size     = 4, color = group_colors[ref_group], fontface = "bold") +
  annotate("text",
           x = max(maaslin_res$coef, na.rm = TRUE),
           y = max(maaslin_res$neg_log10_q, na.rm = TRUE),
           label    = sprintf("n = %d", n_enr_comp),
           hjust    = 1, vjust = 1,
           size     = 4, color = group_colors[comparison_group], fontface = "bold") +
  scale_color_manual(values = volcano_colors) +
  scale_size_continuous(range = c(1.5, 4.5), guide = "none") +
  labs(
    title    = "Differential Abundance — MaAsLin2",
    subtitle = sprintf(
      "%s vs %s  ·  TSS + LOG  ·  BH correction  ·  dashed lines: q = %.2f, |coef| = %.1f",
      comparison_group, ref_group, q_threshold, coef_threshold
    ),
    x     = sprintf(
      "MaAsLin2 Coefficient (log-fold change)\n\u2190 Higher in %s  |  Higher in %s \u2192",
      ref_group, comparison_group
    ),
    y     = expression(-log[10](q-value)),
    color = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold"),
    plot.subtitle   = element_text(color = "grey40", size = 9.5),
    legend.position = "bottom",
    legend.text     = element_text(size = 10)
  )

print(p_volcano)

# ─────────────────────────────────────────────────────────────────────────────
# 9.6  Effect size bar plot — top differentially abundant species
# ─────────────────────────────────────────────────────────────────────────────

n_top_bar <- 20   # max features per direction to show

top_pos <- maaslin_res %>%
  filter(qval < q_threshold, coef > 0) %>%
  arrange(qval, desc(coef)) %>%
  head(n_top_bar / 2)

top_neg <- maaslin_res %>%
  filter(qval < q_threshold, coef < 0) %>%
  arrange(qval, coef) %>%
  head(n_top_bar / 2)

top_features_bar <- bind_rows(top_neg, top_pos) %>%
  mutate(feature_clean = fct_reorder(feature_clean, coef))

if (nrow(top_features_bar) > 0) {
  
  p_bar_da <- ggplot(top_features_bar,
                     aes(x = coef, y = feature_clean, fill = Direction)) +
    geom_col(alpha = 0.85, width = 0.75) +
    geom_errorbarh(
      aes(xmin = coef - stderr, xmax = coef + stderr),
      height    = 0.3,
      color     = "grey25",
      linewidth = 0.5
    ) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
    scale_fill_manual(values = volcano_colors) +
    labs(
      title    = sprintf("Top %d Differentially Abundant Species", nrow(top_features_bar)),
      subtitle = sprintf(
        "MaAsLin2  ·  q < %.2f  ·  Error bars = ±1 SE  ·  TSS + LOG normalization",
        q_threshold
      ),
      x    = sprintf(
        "MaAsLin2 Coefficient\n\u2190 Higher in %s  |  Higher in %s \u2192",
        ref_group, comparison_group
      ),
      y    = NULL,
      fill = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold"),
      plot.subtitle   = element_text(color = "grey40", size = 9.5),
      legend.position = "bottom",
      axis.text.y     = element_text(size = 9)
    )
  
  print(p_bar_da)
  
  # Combined panel: volcano + bar
  p_da_panel <- p_volcano + p_bar_da +
    plot_layout(widths = c(1.4, 1)) +
    plot_annotation(
      title = "Differential Abundance Analysis — MaAsLin2",
      theme = theme(plot.title = element_text(face = "bold", size = 15))
    )
  
  print(p_da_panel)
  
} else {
  cat("\n⚠ No features passed q < 0.25 — bar plot skipped\n")
  cat("  Suggestions:\n")
  cat("  1. Relax threshold: q_threshold <- 0.30\n")
  cat("  2. Check sample sizes per group\n")
  cat("  3. Review prevalence filter in Part 4.7\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# 9.7  Save results table
# ─────────────────────────────────────────────────────────────────────────────

results_export <- maaslin_res %>%
  select(
    Feature       = feature_clean,
    Coefficient   = coef,
    SE            = stderr,
    p_value       = pval,
    q_value       = qval,
    Direction,
    N             = N,
    N_not_zero    = N.not.0
  ) %>%
  arrange(q_value, desc(abs(Coefficient)))

write_csv(results_export,
          file.path(maaslin_output_dir, "maaslin2_results_clean.csv"))

cat(sprintf("\n✔ Clean results table saved: %s/maaslin2_results_clean.csv\n",
            maaslin_output_dir))

# ─────────────────────────────────────────────────────────────────────────────
# 9.8  Full summary
# ─────────────────────────────────────────────────────────────────────────────

cat("\n══════════════════════════════════════════════════════════════════════════\n")
cat("  DIFFERENTIAL ABUNDANCE — MaAsLin2 SUMMARY\n")
cat("══════════════════════════════════════════════════════════════════════════\n")
cat(sprintf("  Reference group         : %s\n",  ref_group))
cat(sprintf("  Comparison group        : %s\n",  comparison_group))
cat(sprintf("  Normalization           : TSS + LOG\n"))
cat(sprintf("  Multiple test correction: Benjamini-Hochberg (BH)\n"))
cat(sprintf("  Features tested         : %d\n",  nrow(maaslin_res)))
cat(sprintf("  Significant (q < 0.25)  : %d\n",  n_sig_25))
cat(sprintf("  Significant (q < 0.05)  : %d\n",  n_sig_05))
cat(sprintf("  Enriched in %-27s: %d\n", comparison_group, n_enr_comp))
cat(sprintf("  Enriched in %-27s: %d\n", ref_group,        n_enr_ref))
cat(sprintf("  Output directory        : %s/\n",  maaslin_output_dir))

# ─────────────────────────────────────────────────────────────────────────────
# PART 10 — EXPORT: PLOTS (PNG) + FORMATTED STATISTICAL TABLES (XLSX)
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# 10.1  Create output directory structure
# ─────────────────────────────────────────────────────────────────────────────

output_root <- "microbiome_results"

export_dirs <- c(
  file.path(output_root, "plots", "01_QC"),
  file.path(output_root, "plots", "02_PCA"),
  file.path(output_root, "plots", "03_Alpha_Diversity"),
  file.path(output_root, "plots", "04_Beta_Diversity"),
  file.path(output_root, "plots", "05_Differential_Abundance"),
  file.path(output_root, "tables")
)

for (d in export_dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("✔ Output directories created under: %s/\n", output_root))

# ─────────────────────────────────────────────────────────────────────────────
# 10.2  Helper functions
# ─────────────────────────────────────────────────────────────────────────────

# Save ggplot as PNG only
save_gg <- function(plot_obj, filepath, width = 10, height = 7, dpi = 300) {
  ggsave(paste0(filepath, ".png"),
         plot   = plot_obj,
         width  = width,
         height = height,
         dpi    = dpi,
         bg     = "white")
  cat(sprintf("  ✔ %s.png  [%.0f×%.0f in, %d dpi]\n",
              basename(filepath), width, height, dpi))
}

# Save ComplexHeatmap as PNG only (requires base graphics device)
save_ht <- function(ht_obj, filepath, width = 14, height = 11, dpi = 300) {
  png(paste0(filepath, ".png"),
      width  = width,
      height = height,
      units  = "in",
      res    = dpi,
      bg     = "white")
  draw(ht_obj, merge_legend = TRUE)
  dev.off()
  cat(sprintf("  ✔ %s.png  [%.0f×%.0f in, %d dpi]\n",
              basename(filepath), width, height, dpi))
}

# Add a formatted sheet to an openxlsx workbook
add_fmt_sheet <- function(wb, sheet_name, data, hdr_color = "#4472C4") {
  
  addWorksheet(wb, sheetName = sheet_name)
  
  # Header style
  hdr_style <- createStyle(
    fontColour     = "#FFFFFF",
    fgFill         = hdr_color,
    fontName       = "Calibri",
    fontSize       = 11,
    textDecoration = "bold",
    halign         = "center",
    valign         = "center",
    border         = "TopBottomLeftRight",
    borderColour   = "#FFFFFF",
    wrapText       = TRUE
  )
  # Body style
  body_style <- createStyle(
    fontName     = "Calibri",
    fontSize     = 10,
    border       = "TopBottomLeftRight",
    borderColour = "#BFBFBF"
  )
  # Alternating row style
  alt_style <- createStyle(
    fontName     = "Calibri",
    fontSize     = 10,
    fgFill       = "#F2F7FF",
    border       = "TopBottomLeftRight",
    borderColour = "#BFBFBF"
  )
  
  writeData(wb, sheet = sheet_name, x = data, headerStyle = hdr_style)
  
  n_rows <- nrow(data)
  n_cols <- ncol(data)
  
  if (n_rows > 0) {
    addStyle(wb, sheet_name,
             style      = body_style,
             rows       = 2:(n_rows + 1),
             cols       = 1:n_cols,
             gridExpand = TRUE)
    
    even_rows <- seq(3, n_rows + 1, by = 2)
    if (length(even_rows) > 0) {
      addStyle(wb, sheet_name,
               style      = alt_style,
               rows       = even_rows,
               cols       = 1:n_cols,
               gridExpand = TRUE)
    }
  }
  
  setColWidths(wb, sheet_name, cols = 1:n_cols, widths = "auto")
  freezePane(wb, sheet_name, firstRow = TRUE)
  setRowHeights(wb, sheet_name, rows = 1, heights = 22)
  
  invisible(wb)
}

# ─────────────────────────────────────────────────────────────────────────────
# 10.3  Export plots — PNG only
# ─────────────────────────────────────────────────────────────────────────────

cat("\n── Exporting plots (PNG, 300 dpi) ───────────────────────────────────────\n")

# ── Part 4: QC ────────────────────────────────────────────────────────────────
d <- file.path(output_root, "plots", "01_QC")
cat(" [QC]\n")
save_gg(p_host_all,     file.path(d, "host_contamination_all_samples"),  width = 13, height = 6)
save_gg(p_host_nonzero, file.path(d, "host_contamination_above_zero"),   width = 11, height = 6)
save_gg(p_depth,        file.path(d, "sequencing_depth_per_sample"),     width = 13, height = 6)

# ── Parts 5–6: PCA & Biplots ──────────────────────────────────────────────────
d <- file.path(output_root, "plots", "02_PCA")
cat(" [PCA]\n")
save_gg(p_pca_12,       file.path(d, "pca_taxonomic_PC1_PC2"),     width = 9,  height = 7)
save_gg(p_pca_13,       file.path(d, "pca_taxonomic_PC1_PC3"),     width = 9,  height = 7)
save_gg(p_scree,        file.path(d, "pca_taxonomic_scree"),       width = 8,  height = 5)
save_gg(p_biplot_12,    file.path(d, "pca_biplot_PC1_PC2"),        width = 11, height = 8)
save_gg(p_biplot_13,    file.path(d, "pca_biplot_PC1_PC3"),        width = 11, height = 8)
save_gg(p_biplot_panel, file.path(d, "pca_biplot_panel"),          width = 20, height = 8)
save_gg(p_pca_ko_12,    file.path(d, "pca_KO_PC1_PC2"),            width = 9,  height = 7)
save_gg(p_pca_ko_13,    file.path(d, "pca_KO_PC1_PC3"),            width = 9,  height = 7)
save_gg(p_scree_ko,     file.path(d, "pca_KO_scree"),              width = 8,  height = 5)
save_gg(p_compare,      file.path(d, "pca_taxonomic_vs_KO"),       width = 16, height = 7)

# ── Part 7: Alpha Diversity ───────────────────────────────────────────────────
d <- file.path(output_root, "plots", "03_Alpha_Diversity")
cat(" [Alpha Diversity]\n")
save_gg(p_depth,        file.path(d, "rarefaction_depth"),         width = 13, height = 6)
save_gg(p_alpha_all,    file.path(d, "alpha_all_metrics_faceted"), width = 12, height = 10)
save_gg(p_alpha_panel,  file.path(d, "alpha_2x2_panel"),           width = 14, height = 12)
save_gg(p_observed,     file.path(d, "alpha_observed_richness"),   width = 7,  height = 6)
save_gg(p_chao1,        file.path(d, "alpha_chao1"),               width = 7,  height = 6)
save_gg(p_shannon,      file.path(d, "alpha_shannon"),             width = 7,  height = 6)
save_gg(p_simpson,      file.path(d, "alpha_simpson"),             width = 7,  height = 6)

# ── Part 8: Beta Diversity ────────────────────────────────────────────────────
d <- file.path(output_root, "plots", "04_Beta_Diversity")
cat(" [Beta Diversity]\n")

pcoa_panel_width <- ifelse(length(pcoa_plots) == 3, 24, 16)

for (dist_name in names(pcoa_plots)) {
  fn <- gsub("[^a-zA-Z0-9]", "_", dist_name)
  save_gg(pcoa_plots[[dist_name]],
          file.path(d, sprintf("pcoa_%s", fn)), width = 9, height = 7)
}
save_gg(p_pcoa_panel, file.path(d, "pcoa_panel_all_distances"),
        width = pcoa_panel_width, height = 7)

for (dist_name in names(nmds_plots)) {
  fn <- gsub("[^a-zA-Z0-9]", "_", dist_name)
  save_gg(nmds_plots[[dist_name]],
          file.path(d, sprintf("nmds_%s", fn)), width = 9, height = 7)
}
save_gg(p_nmds_panel, file.path(d, "nmds_panel_all_distances"),
        width = pcoa_panel_width, height = 7)

# ComplexHeatmap — base graphics device
save_ht(ht, file.path(d, "beta_heatmap_top_species"), width = 14, height = 12)

# ── Part 9: Differential Abundance ───────────────────────────────────────────
d <- file.path(output_root, "plots", "05_Differential_Abundance")
cat(" [Differential Abundance]\n")
save_gg(p_volcano, file.path(d, "volcano_maaslin2"),  width = 11, height = 8)

if (exists("p_bar_da") && !is.null(p_bar_da)) {
  save_gg(p_bar_da,   file.path(d, "barplot_top_da_species"), width = 10, height = 9)
  save_gg(p_da_panel, file.path(d, "da_combined_panel"),      width = 20, height = 8)
}

cat(sprintf("\n✔ All plots saved to: %s/plots/\n", output_root))

# ─────────────────────────────────────────────────────────────────────────────
# 10.4  Export statistical tables — formatted Excel workbooks
# ─────────────────────────────────────────────────────────────────────────────

cat("\n── Exporting statistical tables ─────────────────────────────────────────\n")

tbl_dir <- file.path(output_root, "tables")

# NMDS stress summary (reused across workbooks)
nmds_stress_df <- tibble(
  Distance = names(nmds_stress),
  Stress   = round(unlist(nmds_stress), 4),
  Quality  = case_when(
    Stress < 0.10 ~ "Excellent (< 0.10)",
    Stress < 0.20 ~ "Good (< 0.20)",
    Stress < 0.30 ~ "Acceptable (< 0.30)",
    TRUE          ~ "Poor — interpret with caution"
  )
)

# ── Workbook 1: QC ────────────────────────────────────────────────────────────
wb_qc <- createWorkbook()
add_fmt_sheet(wb_qc, "Host_Contamination",
              host_contamination_df %>%
                arrange(desc(Host_Pct)) %>%
                mutate(Host_Pct = round(Host_Pct, 4)),
              hdr_color = "#7030A0")

saveWorkbook(wb_qc,
             file.path(tbl_dir, "01_QC_host_contamination.xlsx"),
             overwrite = TRUE)
cat("  ✔ 01_QC_host_contamination.xlsx\n")

# ── Workbook 2: Alpha Diversity ───────────────────────────────────────────────
wb_alpha <- createWorkbook()

add_fmt_sheet(wb_alpha, "Descriptive_Statistics",
              alpha_summary,                                hdr_color = "#2E75B6")
add_fmt_sheet(wb_alpha, "Kruskal_Wallis",
              kw_results,                                   hdr_color = "#2E75B6")
add_fmt_sheet(wb_alpha, "Wilcoxon_Rank_Sum",
              wx_results,                                   hdr_color = "#2E75B6")
add_fmt_sheet(wb_alpha, "Per_Sample_Metrics",
              alpha_div %>%
                select(Sample_ID_FM_Pipeline, Group,
                       Observed, Chao1, Shannon, Simpson) %>%
                arrange(Group, Sample_ID_FM_Pipeline),      hdr_color = "#2E75B6")

saveWorkbook(wb_alpha,
             file.path(tbl_dir, "02_alpha_diversity_statistics.xlsx"),
             overwrite = TRUE)
cat("  ✔ 02_alpha_diversity_statistics.xlsx\n")

# ── Workbook 3: Beta Diversity ────────────────────────────────────────────────
wb_beta <- createWorkbook()

add_fmt_sheet(wb_beta, "PERMANOVA_adonis2",  permanova_results, hdr_color = "#375623")
add_fmt_sheet(wb_beta, "ANOSIM",             anosim_results,    hdr_color = "#375623")
add_fmt_sheet(wb_beta, "Betadisper",         betadisp_results,  hdr_color = "#375623")
add_fmt_sheet(wb_beta, "NMDS_Stress",        nmds_stress_df,    hdr_color = "#375623")
add_fmt_sheet(wb_beta, "PCA_Loadings_Top",   loadings_table,    hdr_color = "#375623")

saveWorkbook(wb_beta,
             file.path(tbl_dir, "03_beta_diversity_statistics.xlsx"),
             overwrite = TRUE)
cat("  ✔ 03_beta_diversity_statistics.xlsx\n")

# ── Workbook 4: Differential Abundance ───────────────────────────────────────
wb_da <- createWorkbook()

da_all <- maaslin_res %>%
  select(Feature     = feature_clean,
         Coefficient = coef,
         SE          = stderr,
         p_value     = pval,
         q_value     = qval,
         Direction,
         N,
         N_not_zero  = N.not.0) %>%
  arrange(q_value, desc(abs(Coefficient)))

da_sig          <- da_all %>% filter(q_value < q_threshold)
da_enr_comp     <- da_all %>% filter(q_value < q_threshold, Coefficient > 0)
da_enr_ref      <- da_all %>% filter(q_value < q_threshold, Coefficient < 0)

add_fmt_sheet(wb_da, "All_Results",   da_all,  hdr_color = "#833C00")
add_fmt_sheet(wb_da, "Significant",   da_sig,  hdr_color = "#833C00")

if (nrow(da_enr_comp) > 0)
  add_fmt_sheet(wb_da,
                substr(paste0("Enriched_", comparison_group), 1, 31),
                da_enr_comp,
                hdr_color = group_colors[comparison_group])

if (nrow(da_enr_ref) > 0)
  add_fmt_sheet(wb_da,
                substr(paste0("Enriched_", ref_group), 1, 31),
                da_enr_ref,
                hdr_color = group_colors[ref_group])

saveWorkbook(wb_da,
             file.path(tbl_dir, "04_differential_abundance_maaslin2.xlsx"),
             overwrite = TRUE)
cat("  ✔ 04_differential_abundance_maaslin2.xlsx\n")

# ── Workbook 5: MASTER — all sections in one file ─────────────────────────────
wb_master <- createWorkbook()

add_fmt_sheet(wb_master, "QC_Host_Contamination",
              host_contamination_df %>%
                arrange(desc(Host_Pct)) %>%
                mutate(Host_Pct = round(Host_Pct, 4)),     hdr_color = "#7030A0")
add_fmt_sheet(wb_master, "Alpha_Descriptive",    alpha_summary,     hdr_color = "#2E75B6")
add_fmt_sheet(wb_master, "Alpha_Kruskal_Wallis", kw_results,        hdr_color = "#2E75B6")
add_fmt_sheet(wb_master, "Alpha_Wilcoxon",       wx_results,        hdr_color = "#2E75B6")
add_fmt_sheet(wb_master, "Beta_PERMANOVA",       permanova_results, hdr_color = "#375623")
add_fmt_sheet(wb_master, "Beta_ANOSIM",          anosim_results,    hdr_color = "#375623")
add_fmt_sheet(wb_master, "Beta_Betadisper",      betadisp_results,  hdr_color = "#375623")
add_fmt_sheet(wb_master, "Beta_NMDS_Stress",     nmds_stress_df,    hdr_color = "#375623")
add_fmt_sheet(wb_master, "PCA_Loadings",         loadings_table,    hdr_color = "#375623")
add_fmt_sheet(wb_master, "DA_All_Results",       da_all,            hdr_color = "#833C00")

if (nrow(da_sig) > 0)
  add_fmt_sheet(wb_master, "DA_Significant",     da_sig,            hdr_color = "#833C00")

saveWorkbook(wb_master,
             file.path(tbl_dir, "MASTER_all_results.xlsx"),
             overwrite = TRUE)
cat("  ✔ MASTER_all_results.xlsx  ← all sections combined\n")

# ─────────────────────────────────────────────────────────────────────────────
# 10.5  Export session info
# ─────────────────────────────────────────────────────────────────────────────

session_path <- file.path(output_root, "session_info.txt")
sink(session_path)
cat("══════════════════════════════════════════════════════════════\n")
cat("  MICROBIOME ANALYSIS — SESSION INFORMATION\n")
cat(sprintf("  Generated: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("══════════════════════════════════════════════════════════════\n\n")
print(sessionInfo())
sink()
cat(sprintf("\n✔ Session info saved: session_info.txt\n"))

# ─────────────────────────────────────────────────────────────────────────────
# 10.6  Final export summary
# ─────────────────────────────────────────────────────────────────────────────

cat("\n══════════════════════════════════════════════════════════════════════════\n")
cat("  EXPORT COMPLETE\n")
cat("══════════════════════════════════════════════════════════════════════════\n\n")
cat(sprintf("  Root : %s/\n\n", output_root))
cat("  plots/\n")
cat("  ├── 01_QC/                       host contamination + depth\n")
cat("  ├── 02_PCA/                      PCA, biplots, scree (taxonomic + KO)\n")
cat("  ├── 03_Alpha_Diversity/          violin/boxplots, 2×2 panel\n")
cat("  ├── 04_Beta_Diversity/           PCoA, NMDS, ComplexHeatmap\n")
cat("  └── 05_Differential_Abundance/  volcano, bar plot, panel\n\n")
cat("  tables/\n")
cat("  ├── 01_QC_host_contamination.xlsx\n")
cat("  ├── 02_alpha_diversity_statistics.xlsx\n")
cat("  ├── 03_beta_diversity_statistics.xlsx\n")
cat("  ├── 04_differential_abundance_maaslin2.xlsx\n")
cat("  ├── MASTER_all_results.xlsx          ← all sections combined\n")
cat("  └── session_info.txt\n\n")
cat("  Format : PNG only · 300 dpi · white background\n")
cat("  Tables : formatted Excel · coloured headers · frozen rows\n")
cat("           alternating row shading · auto column widths\n")
cat("══════════════════════════════════════════════════════════════════════════\n")
