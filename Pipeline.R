#Step -by step
# ============================================================
# MICROBIOME ANALYSIS — COMPLETE REVISED CODE v001
# ============================================================

# ─────────────────────────────────────────────────────────────────────────────
# PART 0 — INSTALLATION
# ─────────────────────────────────────────────────────────────────────────────

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "phyloseq", "microbiome", "DESeq2", "ALDEx2",
  "lefser", "KEGGREST", "clusterProfiler",
  "MicrobiomeStat", "Maaslin2"
))

install.packages(c(
  "tidyverse", "vegan", "ggplot2", "ggrepel",
  "reshape2", "RColorBrewer", "pheatmap",
  "ggpubr", "scales", "patchwork"
))

# ─────────────────────────────────────────────────────────────────────────────
# PART 1 — LIBRARIES
# ─────────────────────────────────────────────────────────────────────────────

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







































