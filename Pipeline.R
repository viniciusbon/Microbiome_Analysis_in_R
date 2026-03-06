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
# PART 8 — BETA DIVERSITY: ORDINATIONS, STATS & HEATMAP
# ─────────────────────────────────────────────────────────────────────────────

library(ape)
library(vegan)
library(ComplexHeatmap)
library(circlize)

cat("\n==========================================================\n")
cat("BETA DIVERSITY ANALYSIS\n")
cat("==========================================================\n")

# 8.1  Prepare Phylogenetic Tree (Required for UniFrac) ────────────────────────

# NOTE: UniFrac requires a phylogenetic tree. Since the input data doesn't have one,
# we generate a random tree here SOLELY so the UniFrac function can execute.
# For real biological conclusions using UniFrac, please import your actual .tre file.
if (is.null(phy_tree(ps_rare, errorIfNULL = FALSE))) {
  cat("\n[WARNING] No phylogenetic tree found. Generating a random tree for UniFrac demonstration.\n")
  set.seed(123)
  mock_tree <- rtree(ntaxa(ps_rare), rooted = TRUE, tip.label = taxa_names(ps_rare))
  ps_beta <- merge_phyloseq(ps_rare, mock_tree)
} else {
  ps_beta <- ps_rare
}

# 8.2  Calculate Distance Matrices ─────────────────────────────────────────────

cat("Calculating distance matrices (Bray-Curtis, Jaccard, Weighted UniFrac)...\n")
dist_bray <- phyloseq::distance(ps_beta, method = "bray")
dist_jacc <- phyloseq::distance(ps_beta, method = "jaccard", binary = TRUE)
dist_wuni <- phyloseq::distance(ps_beta, method = "wunifrac")

# Extract metadata for statistical testing
meta_beta <- as(sample_data(ps_beta), "data.frame")

# 8.3  Statistical Tests (PERMANOVA & ANOSIM) ──────────────────────────────────

run_beta_stats <- function(dist_mat, dist_name) {
  # PERMANOVA (adonis2)
  perm_res <- adonis2(dist_mat ~ Group, data = meta_beta, permutations = 999)
  p_perm   <- perm_res$`Pr(>F)`[1]
  r2_perm  <- perm_res$R2[1]
  
  # ANOSIM
  anosim_res <- anosim(dist_mat, meta_beta$Group, permutations = 999)
  p_anosim   <- anosim_res$signif
  r_anosim   <- anosim_res$statistic
  
  cat(sprintf("\n--- %s Distance ---\n", dist_name))
  cat(sprintf("PERMANOVA : R2 = %.4f | p-value = %.4f\n", r2_perm, p_perm))
  cat(sprintf("ANOSIM    : R  = %.4f | p-value = %.4f\n", r_anosim, p_anosim))
  
  return(list(p_perm = p_perm, r2_perm = r2_perm))
}

cat("\n══ Beta Diversity Statistical Tests ════════════════════════════════════\n")
stats_bray <- run_beta_stats(dist_bray, "Bray-Curtis")
stats_jacc <- run_beta_stats(dist_jacc, "Jaccard")
stats_wuni <- run_beta_stats(dist_wuni, "Weighted UniFrac")

# 8.4  Ordinations (PCoA & NMDS) ───────────────────────────────────────────────

# Function to create ordination plot
plot_ordination_custom <- function(ps_obj, dist_mat, ord_method, dist_name, stats) {
  
  ord <- ordinate(ps_obj, method = ord_method, distance = dist_mat)
  
  p <- plot_ordination(ps_obj, ord, color = "Group") +
    geom_point(size = 3.5, alpha = 0.85) +
    stat_ellipse(type = "t", linetype = "dashed", alpha = 0.7, linewidth = 0.8) +
    scale_color_manual(values = group_colors) +
    labs(
      title = sprintf("%s (%s)", ord_method, dist_name),
      subtitle = sprintf("PERMANOVA: p = %.3f, R² = %.3f", stats$p_perm, stats$r2_perm)
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
  return(p)
}

p_pcoa_bray <- plot_ordination_custom(ps_beta, dist_bray, "PCoA", "Bray-Curtis", stats_bray)
p_nmds_bray <- plot_ordination_custom(ps_beta, dist_bray, "NMDS", "Bray-Curtis", stats_bray)
p_pcoa_wuni <- plot_ordination_custom(ps_beta, dist_wuni, "PCoA", "Weighted UniFrac", stats_wuni)

# Display ordination plots using patchwork
p_beta_panel <- (p_pcoa_bray | p_nmds_bray | p_pcoa_wuni) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
print(p_beta_panel)

# ─────────────────────────────────────────────────────────────────────────────
# 8.5  COMPLEXHEATMAP: Top Abundant Taxa
# ─────────────────────────────────────────────────────────────────────────────

cat("\nGenerating ComplexHeatmap for top 40 taxa...\n")

# Transform counts to relative abundance
ps_rel <- transform_sample_counts(ps_beta, function(x) x / sum(x))

# Get the top 40 most abundant taxa across all samples
top40_taxa <- names(sort(taxa_sums(ps_rel), decreasing = TRUE))[1:40]
ps_top40   <- prune_taxa(top40_taxa, ps_rel)

# Extract OTU table as matrix
heat_mat <- as(otu_table(ps_top40), "matrix")

# Ensure taxa are rows and samples are columns
if (!taxa_are_rows(ps_top40)) {
  heat_mat <- t(heat_mat)
}

# Simplify rownames (just keep species name instead of the whole string)
clean_rownames <- sapply(strsplit(rownames(heat_mat), "\\|"), tail, 1)
rownames(heat_mat) <- clean_rownames

# Apply Z-score transformation row-wise (to highlight variance across samples)
heat_mat_z <- t(scale(t(heat_mat)))

# Ensure column names in matrix match metadata precisely
heat_mat_z <- heat_mat_z[, rownames(meta_beta)]

# Create Column Annotation for Group
col_ha <- HeatmapAnnotation(
  Group = meta_beta$Group,
  col   = list(Group = group_colors),
  annotation_name_side = "left"
)

# Build Heatmap
ht <- Heatmap(
  heat_mat_z,
  name              = "Z-Score\n(Rel. Abund)",
  top_annotation    = col_ha,
  show_row_names    = TRUE,
  show_column_names = FALSE,       # Hide sample names to avoid clutter
  row_names_gp      = gpar(fontsize = 9, fontface = "italic"),
  column_title      = "Top 40 Most Abundant Species",
  clustering_distance_columns = "euclidean",
  clustering_method_columns   = "ward.D2",
  clustering_distance_rows    = "pearson",
  row_dend_width    = unit(15, "mm")
)

# Draw Heatmap
draw(ht, merge_legend = TRUE)

cat("✔ Beta Diversity analysis and plotting complete.\n")
