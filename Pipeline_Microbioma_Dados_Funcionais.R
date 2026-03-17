###############################################################################
#                    ╔═══════════════════════════════════════╗
#                    ║   CONFIGURAÇÃO DO USUÁRIO - INÍCIO   ║
#                    ╚═══════════════════════════════════════╝
###############################################################################

data_dir           <- "C:/Users/dti-/Desktop/Padronizados e prontos para analise/HAV1801"
sample_id_col_name <- "Sample_ID"
group_col_name     <- "group"
selected_groups    <- c("control", "glycodex")
age_col_name       <- "age"
selected_ages      <- NULL
group_colors       <- c("#E31A1C", "#1F78B4")
names(group_colors) <- selected_groups

###############################################################################
#                    ╔═══════════════════════════════════════╗
#                    ║    CONFIGURAÇÃO DO USUÁRIO - FIM      ║
#                    ╚═══════════════════════════════════════╝
###############################################################################

library(phyloseq); library(DESeq2); library(edgeR); library(compositions)
library(zCompositions); library(dplyr); library(tibble); library(readr)
library(stringr); library(ggplot2); library(ggrepel); library(corrplot)
library(VennDiagram); library(readxl); library(vegan); library(ape)
library(ggpubr); library(patchwork); library(RColorBrewer)
library(ComplexHeatmap); library(circlize); library(openxlsx)
library(tidyr); library(Maaslin2)

select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
mutate <- dplyr::mutate; arrange <- dplyr::arrange; summarise <- dplyr::summarise

# ═══════════════════════════════════════════════════════════════════════════════
# FASE 1
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n", rep("=", 60), "\n=== FASE 1: IMPORTACAO ===\n", rep("=", 60), "\n")
cat("Diretorio:", data_dir, "\n\n")

patterns <- c(AMR = "AMR_counts", CAZy = "CAZy_counts", COGs = "COGs_count",
              EC = "EC_count", KO = "KO_counts", VF = "VF_counts",
              TAX = "tax_species_taxonomic_profiles", meta = "meta")

data_storage <- list()
for (key in names(patterns)) {
  ff <- list.files(data_dir, pattern = patterns[key], full.names = TRUE, ignore.case = TRUE)
  if (length(ff) == 0) { warning(paste("Not found:", key)); next }
  if (length(ff) > 1) { message(paste("Multiple for", key, "using", basename(ff[1]))); ff <- ff[1] }
  data_storage[[key]] <- if (grepl("\\.xlsx$", ff, ignore.case = TRUE)) read_excel(ff) else read_delim(ff, delim = "\t", show_col_types = FALSE)
  message(paste("Imported:", key, "-", basename(ff)))
}

df_ko <- data_storage$KO; df_AMR <- data_storage$AMR; df_VF <- data_storage$VF
df_CAZy <- data_storage$CAZy; df_EC <- data_storage$EC; df_COGs <- data_storage$COGs
df_TAX <- data_storage$TAX; df_meta <- data_storage$meta

pid <- c("sample","Sample","sample_id","SampleID","sample_ID","Sample_ID")
for (col in pid) {
  if (col %in% colnames(df_meta) && col != sample_id_col_name) {
    colnames(df_meta)[colnames(df_meta) == col] <- sample_id_col_name
    cat(paste("Renomeada:", col, "->", sample_id_col_name, "\n")); break
  } else if (col %in% colnames(df_meta) && col == sample_id_col_name) {
    cat(paste(sample_id_col_name, "ja existe\n")); break
  }
}

cat("\n=== VERIFICACAO ===\n")
for (nm in c("KO","AMR","VF","CAZy","EC","COGs","TAX","META")) {
  obj <- switch(nm, KO=df_ko, AMR=df_AMR, VF=df_VF, CAZy=df_CAZy, EC=df_EC, COGs=df_COGs, TAX=df_TAX, META=df_meta)
  if (!is.null(obj)) cat(paste(nm,":",nrow(obj),"x",ncol(obj),"\n")) else cat(paste(nm,": NAO ENCONTRADO\n"))
}

# ═══════════════════════════════════════════════════════════════════════════════
# FASE 2
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n", rep("=", 60), "\n=== FASE 2: TRANSPOSICAO ===\n", rep("=", 60), "\n")

transpose_functional_data <- function(data, label) {
  fids <- data[[1]]; snames <- names(data)[-1]; cm <- data[, -1]
  td <- as.data.frame(t(cm)); names(td) <- fids
  td[[sample_id_col_name]] <- snames
  td[, c(sample_id_col_name, setdiff(names(td), sample_id_col_name))]
}

for (nm in c("df_ko","df_CAZy","df_COGs","df_EC")) {
  obj <- get(nm, envir = .GlobalEnv)
  if (!is.null(obj)) {
    tr <- transpose_functional_data(obj, nm)
    assign(nm, tr, envir = .GlobalEnv)
    cat(paste("OK", nm, ":", nrow(tr), "x", ncol(tr), "\n"))
  }
}
df_ko <- get("df_ko",.GlobalEnv); df_CAZy <- get("df_CAZy",.GlobalEnv)
df_COGs <- get("df_COGs",.GlobalEnv); df_EC <- get("df_EC",.GlobalEnv)

# FASE 2A — TAX
cat("\n=== FASE 2A: TAX ===\n")
tax_taxonomy_table <- NULL
if (!is.null(df_TAX)) {
  tcands <- c("k","p","c","o","f","g","s","Kingdom","Phylum","Class","Order","Family","Genus","Species")
  tcf <- intersect(names(df_TAX), tcands)
  cat(paste("  Tax cols:", paste(tcf, collapse=", "), "\n"))
  if (length(tcf) > 0) {
    tax_taxonomy_table <- df_TAX[, tcf, drop=FALSE]
    if ("g" %in% tcf && "s" %in% tcf) sl <- paste(df_TAX[["g"]], df_TAX[["s"]], sep="|")
    else if ("s" %in% tcf) sl <- df_TAX[["s"]] else sl <- df_TAX[[tail(tcf,1)]]
    sl <- make.unique(as.character(sl))
    cc <- setdiff(names(df_TAX), tcf)
    td <- as.data.frame(t(as.matrix(df_TAX[, cc, drop=FALSE])))
    colnames(td) <- sl; td[[sample_id_col_name]] <- rownames(td)
    td <- td[, c(sample_id_col_name, setdiff(names(td), sample_id_col_name))]
    rownames(td) <- NULL; df_TAX <- td
    cat(paste("  OK TAX:", nrow(df_TAX), "x", ncol(df_TAX), "\n"))
  } else { df_TAX <- transpose_functional_data(df_TAX, "tax") }
} else { cat("  TAX nao encontrado\n") }

# FASE 2B — PIVOT
cat("\n=== FASE 2B: PIVOT ===\n")
pivot_long_to_wide <- function(dataset, dataset_name, id_col) {
  if (is.null(dataset)) return(NULL)
  if (id_col %in% names(dataset)) {
    nr <- nrow(dataset); nu <- length(unique(dataset[[id_col]]))
    if (nu == nr) { cat(paste("  ", dataset_name, ": largo OK\n")); return(dataset) }
    cat(paste("  ", dataset_name, ": LONGO (", nr, "linhas,", nu, "IDs)\n"))
  } else { return(dataset) }
  pmeta <- c(id_col,"group","age","treatment","study","tissue","age_group","sample")
  nm <- setdiff(names(dataset), pmeta)
  fc <- NULL; for (cl in nm) { if (is.character(dataset[[cl]])||is.factor(dataset[[cl]])) { fc <- cl; break } }
  if (is.null(fc)) fc <- nm[1]
  cat(paste("    Feature col:", fc, "| Unique:", length(unique(dataset[[fc]])), "\n"))
  numcols <- nm[sapply(dataset[, nm, drop=FALSE], is.numeric)]
  ccands <- c("Depth","depth","Hits","hits","Count","count","Score","score")
  vc <- NULL; for (cd in ccands) { if (cd %in% numcols) { vc <- cd; break } }
  if (is.null(vc)) {
    cat("    Contagem de ocorrencias\n")
    pd <- dataset[, c(id_col, fc), drop=FALSE]; pd$count_val <- 1
    ag <- aggregate(as.formula(paste("count_val ~", id_col, "+", fc)), data=pd, FUN=sum)
    wd <- reshape(ag, idvar=id_col, timevar=fc, direction="wide", v.names="count_val")
    names(wd) <- gsub("^count_val\\.","", names(wd))
  } else {
    cat(paste("    Value col:", vc, "\n"))
    pd <- dataset[, c(id_col, fc, vc), drop=FALSE]
    ag <- aggregate(as.formula(paste(vc,"~",id_col,"+",fc)), data=pd, FUN=sum)
    wd <- reshape(ag, idvar=id_col, timevar=fc, direction="wide", v.names=vc)
    pf <- paste0(vc,"."); names(wd) <- gsub(paste0("^",gsub("([.])","\\\\\\1",pf)),"",names(wd))
  }
  wd[is.na(wd)] <- 0; attr(wd,"reshapeWide") <- NULL; rownames(wd) <- NULL
  cat(paste("    Resultado:", nrow(wd), "x", ncol(wd), "\n")); wd
}

df_AMR <- pivot_long_to_wide(df_AMR,"AMR",sample_id_col_name)
df_VF <- pivot_long_to_wide(df_VF,"VF",sample_id_col_name)
df_ko <- pivot_long_to_wide(df_ko,"KO",sample_id_col_name)
df_CAZy <- pivot_long_to_wide(df_CAZy,"CAZy",sample_id_col_name)
df_EC <- pivot_long_to_wide(df_EC,"EC",sample_id_col_name)
df_COGs <- pivot_long_to_wide(df_COGs,"COGs",sample_id_col_name)
df_TAX <- pivot_long_to_wide(df_TAX,"TAX",sample_id_col_name)

cat("=== VERIFICACAO POS-PIVOT ===\n")
for (nm in c("KO","AMR","VF","CAZy","EC","COGs","TAX")) {
  obj <- switch(nm, KO=df_ko, AMR=df_AMR, VF=df_VF, CAZy=df_CAZy, EC=df_EC, COGs=df_COGs, TAX=df_TAX)
  if (!is.null(obj)) { nr <- nrow(obj); nu <- length(unique(obj[[sample_id_col_name]]))
  cat(paste(ifelse(nr==nu,"OK","ERRO"), nm,":", nr,"linhas |",nu,"IDs |",ncol(obj),"cols\n"))
  } else cat(paste("AVISO:",nm,"NULL\n"))
}

# ═══════════════════════════════════════════════════════════════════════════════
# FASE 3
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== FASE 3: PADRONIZACAO NOMES ===\n")
stdnames <- function(x) x %>% str_replace_all("\\.","-") %>% str_replace("_$","")
if (!is.null(df_meta) && sample_id_col_name %in% names(df_meta)) {
  cat("Antes:", paste(head(df_meta[[sample_id_col_name]],3), collapse=", "), "\n")
  df_meta[[sample_id_col_name]] <- stdnames(df_meta[[sample_id_col_name]])
  cat("Depois:", paste(head(df_meta[[sample_id_col_name]],3), collapse=", "), "\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# FASE 4
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== FASE 4: FILTRAGEM PELO METADATA ===\n")
ref_ids <- df_meta[[sample_id_col_name]]
cat("Amostras ref:", length(ref_ids), "\n")
filt <- function(ds, nm, ri) {
  if (is.null(ds)||!sample_id_col_name %in% names(ds)) { cat(paste("AVISO:",nm,"\n")); return(NULL) }
  n0 <- nrow(ds); fd <- ds[ds[[sample_id_col_name]] %in% ri, ]
  cat(paste(" ",nm,":",n0,"->",nrow(fd),"\n")); fd
}
df_ko <- filt(df_ko,"KO",ref_ids); df_AMR <- filt(df_AMR,"AMR",ref_ids)
df_VF <- filt(df_VF,"VF",ref_ids); df_CAZy <- filt(df_CAZy,"CAZy",ref_ids)
df_EC <- filt(df_EC,"EC",ref_ids); df_COGs <- filt(df_COGs,"COGs",ref_ids)
df_TAX <- filt(df_TAX,"TAX",ref_ids)

# ═══════════════════════════════════════════════════════════════════════════════
# FASE 5
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== FASE 5: JOIN METADATA ===\n")
req_cols <- c(sample_id_col_name, group_col_name, age_col_name)
miss <- setdiff(req_cols, names(df_meta))
if (length(miss)>0) stop(paste("ERRO: Colunas faltando:", paste(miss, collapse=", "),
                               "\nDisponiveis:", paste(names(df_meta), collapse=", ")))
msub <- df_meta[, req_cols, drop=FALSE]
cat("Grupos:", paste(unique(df_meta[[group_col_name]]), collapse=", "), "\n")

dfs_list <- list(df_ko=df_ko, df_AMR=df_AMR, df_VF=df_VF, df_CAZy=df_CAZy,
                 df_EC=df_EC, df_COGs=df_COGs, df_TAX=df_TAX)
dfs_j <- lapply(names(dfs_list), function(nm) {
  d <- dfs_list[[nm]]; if (is.null(d)) return(NULL)
  m <- merge(d, msub, by=sample_id_col_name, all.x=TRUE)
  co <- c(sample_id_col_name, group_col_name, age_col_name,
          setdiff(names(m), c(sample_id_col_name, group_col_name, age_col_name)))
  m <- m[, co]; cat(paste("OK",nm,":",ncol(d),"->",ncol(m),"cols\n")); m
})
names(dfs_j) <- names(dfs_list); list2env(dfs_j, envir=.GlobalEnv)
df_ko<-dfs_j$df_ko; df_AMR<-dfs_j$df_AMR; df_VF<-dfs_j$df_VF; df_CAZy<-dfs_j$df_CAZy
df_EC<-dfs_j$df_EC; df_COGs<-dfs_j$df_COGs; df_TAX<-dfs_j$df_TAX

# ═══════════════════════════════════════════════════════════════════════════════
# FASE 6
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== FASE 6: FILTRO GRUPOS ===\n")
cat("Grupos:", paste(selected_groups, collapse=", "), "\n")
dfs_f <- list(df_ko=df_ko, df_AMR=df_AMR, df_VF=df_VF, df_CAZy=df_CAZy,
              df_EC=df_EC, df_COGs=df_COGs, df_TAX=df_TAX)
dfs_g <- lapply(names(dfs_f), function(nm) {
  d <- dfs_f[[nm]]; if (is.null(d)||!group_col_name %in% names(d)) return(NULL)
  fd <- d[d[[group_col_name]] %in% selected_groups, ]
  if (!is.null(selected_ages)) fd <- fd[fd[[age_col_name]] %in% selected_ages, ]
  gt <- table(fd[[group_col_name]])
  cat(paste(" ",nm,":",nrow(d),"->",nrow(fd),"|",paste(names(gt),gt,sep="=",collapse=", "),"\n")); fd
})
names(dfs_g) <- paste0(names(dfs_f),"_grouped"); list2env(dfs_g, envir=.GlobalEnv)
df_ko_grouped<-dfs_g$df_ko_grouped; df_AMR_grouped<-dfs_g$df_AMR_grouped
df_VF_grouped<-dfs_g$df_VF_grouped; df_CAZy_grouped<-dfs_g$df_CAZy_grouped
df_EC_grouped<-dfs_g$df_EC_grouped; df_COGs_grouped<-dfs_g$df_COGs_grouped
df_TAX_grouped<-dfs_g$df_TAX_grouped

# ═══════════════════════════════════════════════════════════════════════════════
# FASE 7 — PCA
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== FASE 7: PCA ===\n")
meta_cols_internal <- c(sample_id_col_name, group_col_name, age_col_name)

do_pca <- function(dataset, dsn) {
  if (is.null(dataset)||!group_col_name %in% names(dataset)) { cat(paste("AVISO:",dsn,"\n")); return(NULL) }
  cat(paste("PCA:",dsn,"...\n"))
  mc <- intersect(meta_cols_internal, names(dataset)); fc <- setdiff(names(dataset), meta_cols_internal)
  ms <- dataset[, mc, drop=FALSE]; fd <- dataset[, fc, drop=FALSE]
  fd <- as.data.frame(lapply(fd, function(x) as.numeric(as.character(x)))); fd[is.na(fd)] <- 0
  vv <- apply(fd,2,var,na.rm=TRUE); vf <- !is.na(vv) & vv>0; fd <- fd[, vf, drop=FALSE]
  cat(paste("  Features:", sum(vf), "\n"))
  if (ncol(fd)<3) { cat("  Poucas features\n"); return(NULL) }
  clr <- compositions::clr(fd+1e-6)
  pr <- prcomp(clr, center=TRUE, scale.=TRUE)
  np <- min(3, ncol(pr$x)); ve <- summary(pr)$importance["Proportion of Variance",1:np]*100
  pdf <- as.data.frame(pr$x[,1:np,drop=FALSE])
  pdf$Sample_ID <- ms[[sample_id_col_name]]; pdf$Group <- ms[[group_col_name]]
  p12 <- ggplot(pdf, aes(PC1,PC2,color=Group,label=Sample_ID)) +
    geom_point(size=3.5,alpha=.85) + geom_text_repel(size=2.5,max.overlaps=30,show.legend=FALSE) +
    stat_ellipse(aes(group=Group),type="norm",linetype="dashed",alpha=.5) +
    scale_color_manual(values=group_colors) +
    labs(title=paste("PCA -",dsn,"(PC1vPC2)"), x=sprintf("PC1(%.1f%%)",ve[1]), y=sprintf("PC2(%.1f%%)",ve[2])) +
    theme_minimal(base_size=13) + theme(plot.title=element_text(face="bold"))
  p13 <- NULL
  if (np>=3) { p13 <- ggplot(pdf, aes(PC1,PC3,color=Group,label=Sample_ID)) +
    geom_point(size=3.5,alpha=.85) + geom_text_repel(size=2.5,max.overlaps=30,show.legend=FALSE) +
    stat_ellipse(aes(group=Group),type="norm",linetype="dashed",alpha=.5) +
    scale_color_manual(values=group_colors) +
    labs(title=paste("PCA -",dsn,"(PC1vPC3)"), x=sprintf("PC1(%.1f%%)",ve[1]), y=sprintf("PC3(%.1f%%)",ve[3])) +
    theme_minimal(base_size=13) + theme(plot.title=element_text(face="bold")) }
  cat(paste("  OK:", paste(sprintf("%.1f%%",ve),collapse=", "), "\n\n"))
  list(pca_result=pr, pca_data=pdf, variance_explained=ve, plot_pc12=p12, plot_pc13=p13,
       dataset_name=dsn, group_distribution=table(pdf$Group))
}

dpca <- list(KO=df_ko_grouped, AMR=df_AMR_grouped, VF=df_VF_grouped,
             CAZy=df_CAZy_grouped, EC=df_EC_grouped, COGs=df_COGs_grouped, TAX=df_TAX_grouped)
pca_results_grouped <- list()
for (dn in names(dpca)) pca_results_grouped[[dn]] <- do_pca(dpca[[dn]], dn)
for (dn in names(pca_results_grouped)) if (!is.null(pca_results_grouped[[dn]])) print(pca_results_grouped[[dn]]$plot_pc12)

variance_summary <- do.call(rbind, lapply(names(pca_results_grouped), function(nm) {
  r <- pca_results_grouped[[nm]]; if (is.null(r)) return(NULL)
  ve <- r$variance_explained; np <- length(ve)
  data.frame(Dataset=nm, PC1=round(ve[1],1), PC2=if(np>=2) round(ve[2],1) else NA,
             PC3=if(np>=3) round(ve[3],1) else NA, Total=round(sum(ve),1), stringsAsFactors=FALSE)
}))
for (g in selected_groups) variance_summary[[paste0(g,"_n")]] <- sapply(variance_summary$Dataset, function(d) {
  gd <- pca_results_grouped[[d]]$group_distribution; if (g %in% names(gd)) as.integer(gd[g]) else 0L })
print(variance_summary)

# ═══════════════════════════════════════════════════════════════════════════════
# FASE 8 — ALPHA DIVERSITY
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n", rep("=",70), "\n=== FASE 8: ALPHA DIVERSITY ===\n", rep("=",70), "\n\n")

datasets_grouped_list <- list(KO=df_ko_grouped, AMR=df_AMR_grouped, VF=df_VF_grouped,
                              CAZy=df_CAZy_grouped, EC=df_EC_grouped, COGs=df_COGs_grouped, TAX=df_TAX_grouped)
datasets_grouped_list <- datasets_grouped_list[!sapply(datasets_grouped_list, is.null)]

output_root <- file.path(data_dir, "functional_diversity_results")
dir.create(output_root, recursive=TRUE, showWarnings=FALSE)

ng <- length(selected_groups)
if (ng<=2) group_palette <- setNames(c("#E31A1C","#1F78B4")[1:ng], selected_groups) else
  group_palette <- setNames(brewer.pal(min(ng,9),"Set1")[1:ng], selected_groups)
if (ng==2) comparisons_list <- list(selected_groups) else comparisons_list <- combn(selected_groups,2,simplify=FALSE)

ecm <- function(ds) {
  fc <- setdiff(names(ds), meta_cols_internal)
  m <- as.data.frame(lapply(ds[,fc], function(x) as.numeric(as.character(x))))
  m[is.na(m)] <- 0; rownames(m) <- ds[[sample_id_col_name]]; as.matrix(m)
}
egv <- function(ds) { g <- ds[[group_col_name]]; names(g) <- ds[[sample_id_col_name]]; g }
ss <- function(p) ifelse(p<.001,"***",ifelse(p<.01,"**",ifelse(p<.05,"*","ns")))
sgg <- function(po, fp, w=10, h=7, d=300) {
  ggsave(paste0(fp,".png"), plot=po, width=w, height=h, dpi=d, bg="white")
  cat(paste("    Salvo:", basename(fp), ".png\n"))
}

alpha_results_all <- list(); alpha_stats_all <- list()
alpha_plots_all <- list(); rarefaction_info_all <- list()

for (ds in names(datasets_grouped_list)) {
  cat(rep("-",60),"\n"); cat(paste("ALPHA:",ds,"\n")); cat(rep("-",60),"\n")
  dat <- datasets_grouped_list[[ds]]; cm <- ecm(dat); gv <- egv(dat)
  cm <- round(cm); cm[cm<0] <- 0; cm <- cm[,colSums(cm)>0,drop=FALSE]
  sd2 <- rowSums(cm); rd <- min(sd2)
  cat(paste("  Matriz:",nrow(cm),"x",ncol(cm),"| Depth:",rd,"\n"))
  rarefaction_info_all[[ds]] <- data.frame(Dataset=ds, Min=rd, Max=max(sd2), Med=median(sd2),
                                           Samples=nrow(cm), Features=ncol(cm), stringsAsFactors=FALSE)
  ddf <- data.frame(Sample=names(sd2), Depth=as.numeric(sd2), Group=gv[names(sd2)], stringsAsFactors=FALSE)
  pd <- ggplot(ddf, aes(reorder(Sample,Depth),Depth,fill=Group)) + geom_col(alpha=.85) +
    geom_hline(yintercept=rd,color="darkred",linetype="dashed") +
    scale_fill_manual(values=group_palette) + labs(title=paste(ds,"Depth"),x="Sample",y="Counts") +
    theme_minimal() + theme(axis.text.x=element_text(angle=90,hjust=1,size=7))
  
  if (rd<10) { cr <- cm } else { set.seed(42); cr <- vegan::rrarefy(cm,rd); cr <- cr[,colSums(cr)>0,drop=FALSE] }
  or <- rowSums(cr>0)
  ch <- tryCatch({t(vegan::estimateR(cr))[,"S.chao1"]}, error=function(e) or)
  sh <- vegan::diversity(cr,"shannon"); si <- vegan::diversity(cr,"simpson")
  adf <- data.frame(Sample=rownames(cr), Group=gv[rownames(cr)], Observed=as.numeric(or),
                    Chao1=as.numeric(ch), Shannon=as.numeric(sh), Simpson=as.numeric(si), stringsAsFactors=FALSE)
  alpha_results_all[[ds]] <- adf
  al <- tidyr::pivot_longer(adf, c("Observed","Chao1","Shannon","Simpson"), names_to="Metric", values_to="Value")
  al$Metric <- factor(al$Metric, levels=c("Observed","Chao1","Shannon","Simpson"))
  
  paf <- ggplot(al, aes(Group,Value,fill=Group,color=Group)) +
    geom_violin(alpha=.35,trim=FALSE) + geom_boxplot(width=.18,alpha=.8,outlier.shape=NA,color="black") +
    geom_jitter(width=.08,size=1.8,alpha=.7,show.legend=FALSE) +
    stat_compare_means(method="wilcox.test",comparisons=comparisons_list,label="p.format",tip.length=.01,size=3.5) +
    scale_fill_manual(values=group_palette) + scale_color_manual(values=group_palette) +
    facet_wrap(~Metric,scales="free_y",nrow=2) +
    labs(title=paste("Alpha -",ds), subtitle=sprintf("Rarefied %d|Wilcoxon",rd), x=NULL, y="Value") +
    theme_bw(base_size=13) + theme(plot.title=element_text(face="bold"),strip.text=element_text(face="bold"),
                                   legend.position="bottom",axis.text.x=element_text(angle=20,hjust=1))
  
  pam <- function(m,yl) { d2<-al[al$Metric==m,]; ggplot(d2,aes(Group,Value,fill=Group,color=Group)) +
    geom_violin(alpha=.3,trim=FALSE) + geom_boxplot(width=.2,alpha=.8,outlier.shape=NA,color="black") +
    geom_jitter(width=.08,size=2.2,alpha=.75,show.legend=FALSE) +
    stat_compare_means(method="wilcox.test",comparisons=comparisons_list,label="p.format",size=3.8) +
    scale_fill_manual(values=group_palette) + scale_color_manual(values=group_palette) +
    labs(title=m,x=NULL,y=yl) + theme_bw(base_size=13) + theme(plot.title=element_text(face="bold",hjust=.5),legend.position="none") }
  po<-pam("Observed","Richness"); pc<-pam("Chao1","Chao1"); ps<-pam("Shannon","Shannon"); pi2<-pam("Simpson","Simpson")
  pap <- (po|pc)/(ps|pi2) + plot_annotation(title=paste("Alpha -",ds))
  
  ml <- c("Observed","Chao1","Shannon","Simpson")
  kw <- do.call(rbind, lapply(ml, function(m) {
    t2<-kruskal.test(adf[[m]]~factor(adf$Group))
    data.frame(Dataset=ds,Metric=m,Statistic=round(t2$statistic,4),p_value=round(t2$p.value,4),Sig=ss(t2$p.value),stringsAsFactors=FALSE)
  }))
  wx <- do.call(rbind, lapply(ml, function(m) { do.call(rbind, lapply(comparisons_list, function(pr) {
    va<-adf[adf$Group==pr[1],m]; vb<-adf[adf$Group==pr[2],m]; t2<-wilcox.test(va,vb,exact=FALSE)
    data.frame(Dataset=ds,Metric=m,Comparison=paste(pr,collapse=" vs "),W=round(t2$statistic,4),
               p_value=round(t2$p.value,4),Sig=ss(t2$p.value),stringsAsFactors=FALSE)
  })) }))
  de <- do.call(rbind, lapply(ml, function(m) { do.call(rbind, lapply(selected_groups, function(g) {
    v<-adf[adf$Group==g,m]
    data.frame(Dataset=ds,Group=g,Metric=m,N=length(v),Mean=round(mean(v),4),Median=round(median(v),4),
               SD=round(sd(v),4),Min=round(min(v),4),Max=round(max(v),4),stringsAsFactors=FALSE)
  })) }))
  
  alpha_stats_all[[ds]] <- list(kruskal_wallis=kw, wilcoxon=wx, descriptive=de)
  alpha_plots_all[[ds]] <- list(depth=pd, facet=paf, panel=pap, observed=po, chao1=pc, shannon=ps, simpson=pi2)
  cat("\n  KW:\n"); print(kw[,c("Metric","Statistic","p_value","Sig")])
  print(pd); print(paf); print(pap)
  cat(paste("\n  OK:",ds,"\n\n"))
}

all_kw <- do.call(rbind, lapply(alpha_stats_all, function(x) x$kruskal_wallis))
all_wx <- do.call(rbind, lapply(alpha_stats_all, function(x) x$wilcoxon))
all_desc <- do.call(rbind, lapply(alpha_stats_all, function(x) x$descriptive))

# ═══════════════════════════════════════════════════════════════════════════════
# FASE 9 — BETA DIVERSITY
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n", rep("=",70), "\n=== FASE 9: BETA DIVERSITY ===\n", rep("=",70), "\n\n")

beta_stats_all <- list(); beta_plots_all <- list()

for (ds in names(datasets_grouped_list)) {
  cat(rep("-",60),"\n"); cat(paste("BETA:",ds,"\n")); cat(rep("-",60),"\n")
  dat <- datasets_grouped_list[[ds]]; cm <- ecm(dat); gv <- egv(dat)
  cm <- round(cm); cm[cm<0]<-0; cm <- cm[,colSums(cm)>0,drop=FALSE]
  rd <- min(rowSums(cm))
  if (rd>=10) { set.seed(42); cr<-vegan::rrarefy(cm,rd); cr<-cr[,colSums(cr)>0,drop=FALSE] } else cr<-cm
  rm2 <- sweep(cr,1,rowSums(cr),"/"); rm2[is.nan(rm2)]<-0
  mb <- data.frame(Sample=rownames(cr), Group=gv[rownames(cr)], stringsAsFactors=FALSE)
  cat(paste("  Matriz:",nrow(cr),"x",ncol(cr),"| Rare:",rd,"\n"))
  db <- vegdist(rm2,"bray"); dj <- vegdist(cr,"jaccard",binary=TRUE)
  dl <- list("Bray-Curtis"=db, "Jaccard"=dj)
  
  sbd <- list()
  for (dn in names(dl)) {
    set.seed(42); pm<-adonis2(dl[[dn]]~Group,data=mb,permutations=999,by="margin")
    pR2<-round(pm$R2[1],3); pF2<-round(pm$F[1],3); pp<-pm$`Pr(>F)`[1]
    set.seed(42); an<-anosim(dl[[dn]],mb$Group,permutations=999); aR<-round(an$statistic,3); ap<-an$signif
    at<-sprintf("PERMANOVA: R2=%.3f p=%.4f %s\nANOSIM: R=%.3f p=%.4f %s",pR2,pp,ss(pp),aR,ap,ss(ap))
    cat(sprintf("  [%s] R2=%.3f p=%.4f%s | R=%.3f p=%.4f%s\n",dn,pR2,pp,ss(pp),aR,ap,ss(ap)))
    sbd[[dn]] <- list(permanova_R2=pR2,permanova_F=pF2,permanova_p=pp,permanova_sig=ss(pp),
                      anosim_R=aR,anosim_p=ap,anosim_sig=ss(ap),annot_text=at)
  }
  
  bdr <- do.call(rbind, lapply(names(dl), function(dn) {
    bd<-betadisper(dl[[dn]],mb$Group); set.seed(42); pd2<-permutest(bd,permutations=999)
    data.frame(Dataset=ds,Distance=dn,F_stat=round(pd2$tab[1,"F"],4),
               p_value=round(pd2$tab[1,"Pr(>F)"],4),stringsAsFactors=FALSE) }))
  
  pp_list <- list()
  for (dn in names(dl)) {
    pr2<-pcoa(dl[[dn]]); ei<-pr2$values$Eigenvalues; vp<-round(ei[ei>0]/sum(ei[ei>0])*100,1)
    pd2<-as.data.frame(pr2$vectors[,1:2]); colnames(pd2)<-c("A1","A2")
    pd2$Sample<-rownames(pd2); pd2$Group<-mb$Group[match(pd2$Sample,mb$Sample)]
    p<-ggplot(pd2,aes(A1,A2,color=Group,label=Sample))+geom_point(size=3.5,alpha=.85)+
      geom_text_repel(size=2.5,max.overlaps=25,show.legend=FALSE)+
      stat_ellipse(aes(group=Group),type="norm",linetype="dashed",alpha=.6)+
      scale_color_manual(values=group_palette)+
      annotate("label",x=Inf,y=-Inf,label=sbd[[dn]]$annot_text,hjust=1.03,vjust=-.15,size=2.8,fill="white",alpha=.88,family="mono")+
      labs(title=sprintf("PCoA-%s-%s",ds,dn),x=sprintf("Axis1(%.1f%%)",vp[1]),y=sprintf("Axis2(%.1f%%)",vp[2]))+
      theme_bw(base_size=13)+theme(plot.title=element_text(face="bold"))
    pp_list[[dn]]<-p; print(p)
  }
  pcp<-wrap_plots(pp_list,ncol=2)+plot_layout(guides="collect")+
    plot_annotation(title=paste("PCoA -",ds))&theme(legend.position="bottom"); print(pcp)
  
  np_list<-list(); nst<-list()
  for (dn in names(dl)) {
    set.seed(42); nr2<-tryCatch(metaMDS(dl[[dn]],k=2,trymax=100,trace=FALSE,autotransform=FALSE),error=function(e) NULL)
    if (is.null(nr2)) next; sv<-nr2$stress; nst[[dn]]<-sv
    sl2<-ifelse(sv<.1,"excelente",ifelse(sv<.2,"bom",ifelse(sv<.3,"aceitavel","ruim")))
    ns2<-as.data.frame(scores(nr2,display="sites")); ns2$Sample<-rownames(ns2)
    ns2$Group<-mb$Group[match(ns2$Sample,mb$Sample)]
    p<-ggplot(ns2,aes(NMDS1,NMDS2,color=Group,label=Sample))+geom_point(size=3.5,alpha=.85)+
      geom_text_repel(size=2.5,max.overlaps=25,show.legend=FALSE)+
      stat_ellipse(aes(group=Group),type="norm",linetype="dashed",alpha=.6)+
      scale_color_manual(values=group_palette)+
      annotate("text",x=-Inf,y=Inf,label=sprintf("Stress:%.4f(%s)",sv,sl2),hjust=-.1,vjust=1.4,size=3.5,color="grey30",fontface="italic")+
      labs(title=sprintf("NMDS-%s-%s",ds,dn),x="NMDS1",y="NMDS2")+
      theme_bw(base_size=13)+theme(plot.title=element_text(face="bold"))
    np_list[[dn]]<-p; print(p)
  }
  pnp<-NULL; if(length(np_list)>0) { pnp<-wrap_plots(np_list,ncol=2)+plot_layout(guides="collect")+
    plot_annotation(title=paste("NMDS-",ds))&theme(legend.position="bottom"); print(pnp) }
  
  nth<-min(50,ncol(rm2)); fv2<-apply(rm2,2,var); tf2<-names(sort(fv2,decreasing=TRUE))[1:nth]
  hm2<-t(rm2[,tf2]); hz<-t(scale(t(hm2))); hz<-pmax(pmin(hz,3),-3); hz[is.nan(hz)]<-0
  hcc<-tryCatch(hclust(db,"ward.D2"),error=function(e) hclust(dist(t(hz)),"ward.D2"))
  hcr<-hclust(dist(hz),"ward.D2")
  cf2<-colorRamp2(c(-3,-1.5,0,1.5,3),c("#2166AC","#92C5DE","white","#F4A582","#D6604D"))
  ca2<-HeatmapAnnotation(Group=mb$Group,col=list(Group=group_palette))
  rfs<-ifelse(nth>40,5.5,ifelse(nth>25,7,8.5))
  ht2<-Heatmap(hz,name="Z",col=cf2,top_annotation=ca2,cluster_rows=hcr,cluster_columns=hcc,
               show_column_names=TRUE,column_names_gp=gpar(fontsize=7),show_row_names=TRUE,
               row_names_gp=gpar(fontsize=rfs),row_names_side="left",
               column_title=sprintf("%s-Top%d",ds,nth),border=TRUE)
  draw(ht2,merge_legend=TRUE)
  
  pt2<-do.call(rbind,lapply(names(sbd),function(dn){s<-sbd[[dn]]
  data.frame(Dataset=ds,Distance=dn,R2=s$permanova_R2,F_stat=s$permanova_F,
             p_value=s$permanova_p,Sig=s$permanova_sig,stringsAsFactors=FALSE)}))
  at2<-do.call(rbind,lapply(names(sbd),function(dn){s<-sbd[[dn]]
  data.frame(Dataset=ds,Distance=dn,R_stat=s$anosim_R,p_value=s$anosim_p,Sig=s$anosim_sig,stringsAsFactors=FALSE)}))
  nt2<-if(length(nst)>0) data.frame(Dataset=ds,Distance=names(nst),Stress=round(unlist(nst),4),stringsAsFactors=FALSE) else
    data.frame(Dataset=character(),Distance=character(),Stress=numeric(),stringsAsFactors=FALSE)
  
  beta_stats_all[[ds]]<-list(permanova=pt2,anosim=at2,betadisper=bdr,nmds_stress=nt2)
  beta_plots_all[[ds]]<-list(pcoa_plots=pp_list,pcoa_panel=pcp,nmds_plots=np_list,nmds_panel=pnp,heatmap=ht2)
  cat(paste("\n  OK:",ds,"\n\n"))
}

all_permanova<-do.call(rbind,lapply(beta_stats_all,function(x) x$permanova))
all_anosim<-do.call(rbind,lapply(beta_stats_all,function(x) x$anosim))
all_betadisp<-do.call(rbind,lapply(beta_stats_all,function(x) x$betadisper))
all_nmds_str<-do.call(rbind,lapply(beta_stats_all,function(x) x$nmds_stress))
cat("PERMANOVA:\n"); print(all_permanova)
cat("\nANOSIM:\n"); print(all_anosim)

# ═══════════════════════════════════════════════════════════════════════════════
# FASE 11 — DIFFERENTIAL ABUNDANCE: MaAsLin2 + VOLCANO
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n", rep("=",70), "\n=== FASE 11: MAASLIN2 ===\n", rep("=",70), "\n\n")

q_threshold <- 0.25; coef_threshold <- 0.5; n_top_bar <- 20
ref_group <- selected_groups[1]; comparison_group <- selected_groups[2]
cat(sprintf("Ref: %s | Comp: %s | q: %.2f\n\n", ref_group, comparison_group, q_threshold))

volcano_colors <- c(
  setNames(group_palette[comparison_group], paste("Enriched in", comparison_group)),
  setNames(group_palette[ref_group], paste("Enriched in", ref_group)),
  "Not significant" = "grey72")

da_results_all <- list(); da_plots_all <- list(); da_stats_all <- list()

for (ds in names(datasets_grouped_list)) {
  cat(rep("-",60),"\n"); cat(paste("DA:",ds,"\n")); cat(rep("-",60),"\n")
  dat <- datasets_grouped_list[[ds]]
  fc <- setdiff(names(dat), meta_cols_internal)
  mf <- as.data.frame(lapply(dat[,fc], function(x) as.numeric(as.character(x))))
  mf[is.na(mf)] <- 0; rownames(mf) <- dat[[sample_id_col_name]]
  mf <- mf[, colSums(mf)>0, drop=FALSE]
  mm <- data.frame(Group=dat[[group_col_name]], stringsAsFactors=FALSE)
  rownames(mm) <- dat[[sample_id_col_name]]
  mm$Group <- relevel(factor(mm$Group), ref=ref_group)
  cat(sprintf("  Input: %d samples x %d features\n", nrow(mf), ncol(mf)))
  
  mod <- file.path(output_root, ds, "maaslin2_output")
  dir.create(mod, recursive=TRUE, showWarnings=FALSE)
  
  set.seed(42)
  fit <- tryCatch(Maaslin2(input_data=mf, input_metadata=mm, output=mod, fixed_effects="Group",
                           normalization="TSS", transform="LOG", analysis_method="LM",
                           min_abundance=0, min_prevalence=0, max_significance=0.25,
                           correction="BH", standardize=FALSE, plot_heatmap=FALSE,
                           plot_scatter=FALSE, cores=1),
                  error=function(e) { cat(paste("  ERRO:",e$message,"\n")); NULL })
  if (is.null(fit)) { cat("  Pulando.\n\n"); next }
  
  mr <- fit$results; mr <- mr[mr$metadata=="Group",]
  mr <- mr[order(mr$qval, mr$pval),]
  mr$feature_clean <- sub("^[^|]*\\|","", mr$feature)
  mr$qval_safe <- pmax(mr$qval, 1e-10); mr$neg_log10_q <- -log10(mr$qval_safe)
  mr$Direction <- ifelse(mr$qval<q_threshold & mr$coef>0, paste("Enriched in",comparison_group),
                         ifelse(mr$qval<q_threshold & mr$coef<0, paste("Enriched in",ref_group), "Not significant"))
  
  n25 <- sum(mr$qval<.25,na.rm=TRUE); n05 <- sum(mr$qval<.05,na.rm=TRUE)
  nec <- sum(mr$qval<q_threshold & mr$coef>0, na.rm=TRUE)
  ner <- sum(mr$qval<q_threshold & mr$coef<0, na.rm=TRUE)
  cat(sprintf("  Tested:%d | q<.25:%d | q<.05:%d | Enr %s:%d | Enr %s:%d\n",
              nrow(mr),n25,n05,comparison_group,nec,ref_group,ner))
  
  sf <- mr[mr$qval<q_threshold,]
  tl <- do.call(rbind, lapply(unique(sf$Direction), function(d) {
    s2<-sf[sf$Direction==d,]; s2<-s2[order(s2$qval,-abs(s2$coef)),]; head(s2,12) }))
  
  pv <- ggplot(mr, aes(coef, neg_log10_q, color=Direction, size=abs(coef))) +
    geom_point(alpha=.72) +
    geom_hline(yintercept=-log10(q_threshold), linetype="dashed", color="grey35") +
    geom_vline(xintercept=c(-coef_threshold,coef_threshold), linetype="dashed", color="grey35") +
    scale_color_manual(values=volcano_colors) + scale_size_continuous(range=c(1.5,4.5),guide="none") +
    labs(title=paste("DA -",ds,"- MaAsLin2"),
         subtitle=sprintf("%s vs %s | TSS+LOG | BH | q=%.2f",comparison_group,ref_group,q_threshold),
         x=sprintf("Coefficient\n<- %s | %s ->",ref_group,comparison_group),
         y=expression(-log[10](q)), color=NULL) +
    theme_bw(base_size=13) + theme(plot.title=element_text(face="bold"),
                                   plot.subtitle=element_text(color="grey40",size=9.5),
                                   legend.position="bottom")
  
  if (!is.null(tl) && nrow(tl)>0) pv <- pv +
    geom_text_repel(data=tl, aes(label=feature_clean), size=2.9, max.overlaps=30,
                    segment.color="grey50", box.padding=.45, show.legend=FALSE)
  if (nrow(mr)>0) pv <- pv +
    annotate("text",x=min(mr$coef,na.rm=TRUE),y=max(mr$neg_log10_q,na.rm=TRUE),
             label=sprintf("n=%d",ner),hjust=0,vjust=1,size=4,color=group_palette[ref_group],fontface="bold") +
    annotate("text",x=max(mr$coef,na.rm=TRUE),y=max(mr$neg_log10_q,na.rm=TRUE),
             label=sprintf("n=%d",nec),hjust=1,vjust=1,size=4,color=group_palette[comparison_group],fontface="bold") +
    annotate("text",x=max(mr$coef,na.rm=TRUE),y=-log10(q_threshold),
             label=sprintf("q=%.2f",q_threshold),hjust=1,vjust=-.4,size=3.2,color="grey35",fontface="italic")
  print(pv)
  
  pb <- NULL; ppnl <- NULL
  tp <- head(mr[mr$qval<q_threshold & mr$coef>0,], n_top_bar/2)
  tn <- head(mr[mr$qval<q_threshold & mr$coef<0,], n_top_bar/2)
  tb <- rbind(tn, tp)
  if (nrow(tb)>0) {
    tb$feature_clean <- factor(tb$feature_clean, levels=tb$feature_clean[order(tb$coef)])
    pb <- ggplot(tb, aes(coef,feature_clean,fill=Direction)) +
      geom_col(alpha=.85,width=.75) +
      geom_errorbarh(aes(xmin=coef-stderr,xmax=coef+stderr),height=.3,color="grey25") +
      geom_vline(xintercept=0) + scale_fill_manual(values=volcano_colors) +
      labs(title=sprintf("Top %d DA - %s",nrow(tb),ds),
           subtitle=sprintf("q<%.2f | +/-1SE",q_threshold),
           x=sprintf("Coef\n<- %s | %s ->",ref_group,comparison_group), y=NULL, fill=NULL) +
      theme_bw(base_size=12) + theme(plot.title=element_text(face="bold"),legend.position="bottom")
    print(pb)
    ppnl <- pv + pb + plot_layout(widths=c(1.4,1)) +
      plot_annotation(title=paste("DA -",ds)); print(ppnl)
  }
  
  dex <- data.frame(Dataset=ds, Feature=mr$feature_clean, Coefficient=round(mr$coef,4),
                    SE=round(mr$stderr,4), p_value=round(mr$pval,6), q_value=round(mr$qval,6),
                    Direction=mr$Direction, N=mr$N, N_not_zero=mr$N.not.0, stringsAsFactors=FALSE)
  da_results_all[[ds]] <- dex
  da_stats_all[[ds]] <- data.frame(Dataset=ds, Tested=nrow(mr), Sig_q025=n25, Sig_q005=n05,
                                   Enr_comp=nec, Enr_ref=ner, stringsAsFactors=FALSE)
  da_plots_all[[ds]] <- list(volcano=pv, bar=pb, panel=ppnl)
  cat(paste("\n  OK:",ds,"\n\n"))
}

all_da_stats <- do.call(rbind, da_stats_all)
all_da_results <- do.call(rbind, da_results_all)
all_da_sig <- if(!is.null(all_da_results)) all_da_results[all_da_results$q_value<q_threshold,] else data.frame()
cat("\n=== RESUMO DA ===\n"); print(all_da_stats)

library(Maaslin2)   # BiocManager::install("Maaslin2")
library(forcats)


# ═══════════════════════════════════════════════════════════════════════════════
# FASE 11 — ANÁLISE DE ABUNDÂNCIA DIFERENCIAL (MaAsLin2 + Volcano Plots)
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n", rep("=", 70), "\n")
cat("=== FASE 11: ANALISE DE ABUNDANCIA DIFERENCIAL ===\n")
cat(rep("=", 70), "\n\n")

# ─── 11.0  Verificação de bibliotecas ────────────────────────────────────────

if (!requireNamespace("Maaslin2", quietly = TRUE)) {
  stop("Pacote 'Maaslin2' nao encontrado.\n",
       "  Instale com: BiocManager::install('Maaslin2')")
}
library(Maaslin2)
if (!requireNamespace("forcats", quietly = TRUE))
  install.packages("forcats", quiet = TRUE)
library(forcats)

# ─── 11.1  Configuração ──────────────────────────────────────────────────────

abundance_q_threshold    <- 0.25    # limiar FDR (q-value) — padrão MaAsLin2
abundance_p_threshold    <- 0.05    # limiar p-value para volcano alternativo
abundance_coef_threshold <- 0.5     # tamanho mínimo de efeito (linhas tracejadas)
abundance_n_top_labels   <- 15      # rótulos por direção no volcano
abundance_n_top_bar      <- 30      # features máximas no barplot

abundance_ref_group  <- selected_groups[1]   # grupo de referência
abundance_comp_group <- selected_groups[2]   # grupo de comparação

cat(sprintf("  Grupo de referencia  : %s\n",   abundance_ref_group))
cat(sprintf("  Grupo de comparacao  : %s\n",   abundance_comp_group))
cat(sprintf("  Limiar q-value (FDR) : %.2f\n", abundance_q_threshold))
cat(sprintf("  Limiar p-value       : %.2f\n", abundance_p_threshold))
cat(sprintf("  Limiar |coef|        : %.1f\n", abundance_coef_threshold))
cat(sprintf("  Datasets             : %s\n\n",
            paste(names(datasets_grouped_list), collapse = ", ")))

# ─── 11.2  Diretórios de saída ───────────────────────────────────────────────

abundance_output_root <- file.path(output_root, "abundance")
dir.create(abundance_output_root, recursive = TRUE, showWarnings = FALSE)

for (ds_name in names(datasets_grouped_list)) {
  dir.create(file.path(abundance_output_root, ds_name),
             recursive = TRUE, showWarnings = FALSE)
}

# ─── 11.3  Paleta do volcano ─────────────────────────────────────────────────

abundance_volcano_colors <- c(
  setNames(group_colors[abundance_comp_group],
           paste("Enriched in", abundance_comp_group)),
  setNames(group_colors[abundance_ref_group],
           paste("Enriched in", abundance_ref_group)),
  "Not significant" = "grey72"
)

# ─── 11.4  Armazenamento ─────────────────────────────────────────────────────

da_results_all       <- list()
da_plots_all         <- list()
da_summary_all       <- list()
da_heatmaps_all      <- list()
da_heatmap_nfeats    <- list()

# ═══════════════════════════════════════════════════════════════════════════════
# 11.5  LOOP PRINCIPAL — MaAsLin2 + Volcano + BarPlot + Heatmap
# ═══════════════════════════════════════════════════════════════════════════════

for (ds_name in names(datasets_grouped_list)) {
  
  cat(rep("-", 65), "\n")
  cat(sprintf("  ABUNDANCIA DIFERENCIAL: %s\n", ds_name))
  cat(rep("-", 65), "\n")
  
  dataset <- datasets_grouped_list[[ds_name]]
  
  # ── Preparar matriz de features ─────────────────────────────────────────────
  feature_cols <- setdiff(names(dataset), meta_cols_internal)
  features_df  <- dataset[, feature_cols, drop = FALSE]
  features_df  <- as.data.frame(
    lapply(features_df, function(x) as.numeric(as.character(x)))
  )
  features_df[is.na(features_df)] <- 0
  rownames(features_df) <- dataset[[sample_id_col_name]]
  
  # Remover colunas com soma zero ou variância zero
  col_sums   <- colSums(features_df)
  col_vars   <- apply(features_df, 2, var, na.rm = TRUE)
  valid_cols <- (col_sums > 0) & (!is.na(col_vars)) & (col_vars > 0)
  n_removed  <- sum(!valid_cols)
  features_df <- features_df[, valid_cols, drop = FALSE]
  
  cat(sprintf("  Matriz: %d amostras x %d features (removidas: %d)\n",
              nrow(features_df), ncol(features_df), n_removed))
  
  if (ncol(features_df) < 5) {
    cat(sprintf("  AVISO: Poucas features (%d). Pulando %s.\n\n",
                ncol(features_df), ds_name))
    next
  }
  
  # ── Preparar metadata ──────────────────────────────────────────────────────
  meta_df <- data.frame(
    Group = factor(dataset[[group_col_name]],
                   levels = c(abundance_ref_group, abundance_comp_group)),
    row.names = dataset[[sample_id_col_name]],
    stringsAsFactors = FALSE
  )
  
  n_ref  <- sum(meta_df$Group == abundance_ref_group,  na.rm = TRUE)
  n_comp <- sum(meta_df$Group == abundance_comp_group, na.rm = TRUE)
  cat(sprintf("  %s: n=%d  |  %s: n=%d\n",
              abundance_ref_group, n_ref, abundance_comp_group, n_comp))
  
  if (n_ref < 3 || n_comp < 3) {
    cat(sprintf("  AVISO: <3 amostras num grupo. Pulando %s.\n\n", ds_name))
    next
  }
  
  # ── Executar MaAsLin2 ──────────────────────────────────────────────────────
  maaslin_out_dir <- file.path(abundance_output_root, ds_name,
                               "maaslin2_raw_output")
  dir.create(maaslin_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat("  Executando MaAsLin2... ")
  
  set.seed(42)
  fit <- tryCatch(
    Maaslin2(
      input_data       = features_df,
      input_metadata   = meta_df,
      output           = maaslin_out_dir,
      fixed_effects    = "Group",
      normalization    = "TSS",
      transform        = "LOG",
      analysis_method  = "LM",
      min_abundance    = 0,
      min_prevalence   = 0,
      max_significance = 0.25,
      correction       = "BH",
      standardize      = FALSE,
      plot_heatmap     = FALSE,
      plot_scatter     = FALSE,
      cores            = 1
    ),
    error = function(e) {
      cat(sprintf("\n  ERRO MaAsLin2: %s\n", e$message))
      NULL
    }
  )
  
  if (is.null(fit) || is.null(fit$results) || nrow(fit$results) == 0) {
    cat(sprintf("\n  AVISO: Sem resultados para %s. Pulando.\n\n", ds_name))
    next
  }
  cat("OK!\n")
  
  # ── Processar resultados ────────────────────────────────────────────────────
  res <- fit$results %>%
    as_tibble() %>%
    filter(metadata == "Group") %>%
    arrange(qval, pval) %>%
    mutate(
      Dataset       = ds_name,
      feature_clean = feature,
      qval_safe     = pmax(qval, 1e-300),
      pval_safe     = pmax(pval, 1e-300),
      neg_log10_q   = -log10(qval_safe),
      neg_log10_p   = -log10(pval_safe),
      Direction = case_when(
        qval < abundance_q_threshold & coef > 0 ~
          paste("Enriched in", abundance_comp_group),
        qval < abundance_q_threshold & coef < 0 ~
          paste("Enriched in", abundance_ref_group),
        TRUE ~ "Not significant"
      ),
      Direction_pval = case_when(
        pval < abundance_p_threshold & coef > 0 ~
          paste("Enriched in", abundance_comp_group),
        pval < abundance_p_threshold & coef < 0 ~
          paste("Enriched in", abundance_ref_group),
        TRUE ~ "Not significant"
      )
    )
  
  # ── Contagens ───────────────────────────────────────────────────────────────
  n_total   <- nrow(res)
  n_sig_q25 <- sum(res$qval < 0.25, na.rm = TRUE)
  n_sig_q10 <- sum(res$qval < 0.10, na.rm = TRUE)
  n_sig_q05 <- sum(res$qval < 0.05, na.rm = TRUE)
  n_sig_p05 <- sum(res$pval < 0.05, na.rm = TRUE)
  n_up_q    <- sum(res$qval < abundance_q_threshold & res$coef > 0, na.rm = TRUE)
  n_down_q  <- sum(res$qval < abundance_q_threshold & res$coef < 0, na.rm = TRUE)
  n_up_p    <- sum(res$pval < abundance_p_threshold & res$coef > 0, na.rm = TRUE)
  n_down_p  <- sum(res$pval < abundance_p_threshold & res$coef < 0, na.rm = TRUE)
  
  cat(sprintf("\n  Resultados:\n"))
  cat(sprintf("    Features testadas        : %d\n", n_total))
  cat(sprintf("    Significativas (q < 0.25): %d\n", n_sig_q25))
  cat(sprintf("    Significativas (q < 0.10): %d\n", n_sig_q10))
  cat(sprintf("    Significativas (q < 0.05): %d\n", n_sig_q05))
  cat(sprintf("    Significativas (p < 0.05): %d\n", n_sig_p05))
  cat(sprintf("    Enriquecidas %s (q)      : %d\n", abundance_comp_group, n_up_q))
  cat(sprintf("    Enriquecidas %s (q)      : %d\n", abundance_ref_group,  n_down_q))
  
  if (n_sig_q25 > 0) {
    cat("\n    Top features (q < 0.25):\n")
    print(
      res %>%
        filter(qval < 0.25) %>%
        select(feature_clean, coef, stderr, pval, qval, Direction) %>%
        head(20),
      n = 20
    )
  }
  
  # ── Labels para volcano ────────────────────────────────────────────────────
  top_labels_q <- res %>%
    filter(qval < abundance_q_threshold) %>%
    group_by(Direction) %>%
    arrange(qval, desc(abs(coef))) %>%
    slice_head(n = abundance_n_top_labels) %>%
    ungroup()
  
  top_labels_p <- res %>%
    filter(pval < abundance_p_threshold) %>%
    group_by(Direction_pval) %>%
    arrange(pval, desc(abs(coef))) %>%
    slice_head(n = abundance_n_top_labels) %>%
    ungroup()
  
  x_range <- range(res$coef, na.rm = TRUE)
  y_max_q <- max(res$neg_log10_q, na.rm = TRUE)
  y_max_p <- max(res$neg_log10_p, na.rm = TRUE)
  
  # ═════════════════════════════════════════════════════════════════════════════
  #  VOLCANO PLOT  (q-value)
  # ═════════════════════════════════════════════════════════════════════════════
  
  p_volc_q <- ggplot(res, aes(x = coef, y = neg_log10_q,
                              color = Direction, size = abs(coef))) +
    geom_point(alpha = 0.72) +
    geom_hline(yintercept = -log10(abundance_q_threshold),
               linetype = "dashed", color = "grey35", linewidth = 0.6) +
    geom_vline(xintercept = c(-abundance_coef_threshold,
                              abundance_coef_threshold),
               linetype = "dashed", color = "grey35", linewidth = 0.6) +
    scale_color_manual(values = abundance_volcano_colors) +
    scale_size_continuous(range = c(1.5, 4.5), guide = "none") +
    annotate("text",
             x = x_range[2] * 0.95, y = -log10(abundance_q_threshold),
             label    = sprintf("q = %.2f", abundance_q_threshold),
             hjust = 1, vjust = -0.5, size = 3.2,
             color = "grey35", fontface = "italic") +
    annotate("text",
             x = x_range[1], y = y_max_q,
             label    = sprintf("n = %d", n_down_q),
             hjust = 0, vjust = 1, size = 4,
             color = unname(group_colors[abundance_ref_group]),
             fontface = "bold") +
    annotate("text",
             x = x_range[2], y = y_max_q,
             label    = sprintf("n = %d", n_up_q),
             hjust = 1, vjust = 1, size = 4,
             color = unname(group_colors[abundance_comp_group]),
             fontface = "bold") +
    labs(
      title = sprintf("Volcano Plot \u2014 %s (q-value)", ds_name),
      subtitle = sprintf(
        "%s vs %s  \u00B7  MaAsLin2: TSS+LOG  \u00B7  BH correction  \u00B7  q < %.2f",
        abundance_comp_group, abundance_ref_group, abundance_q_threshold),
      x = sprintf("MaAsLin2 Coefficient\n\u2190 Higher in %s  |  Higher in %s \u2192",
                  abundance_ref_group, abundance_comp_group),
      y     = expression(-log[10](q-value)),
      color = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(color = "grey40", size = 9.5),
          legend.position = "bottom",
          legend.text     = element_text(size = 10))
  
  if (nrow(top_labels_q) > 0) {
    p_volc_q <- p_volc_q +
      geom_text_repel(data          = top_labels_q,
                      aes(label     = feature_clean),
                      size          = 2.7,
                      max.overlaps  = 30,
                      segment.color = "grey50",
                      segment.size  = 0.3,
                      box.padding   = 0.45,
                      show.legend   = FALSE)
  }
  
  # ═════════════════════════════════════════════════════════════════════════════
  #  VOLCANO PLOT  (p-value)
  # ═════════════════════════════════════════════════════════════════════════════
  
  p_volc_p <- ggplot(res, aes(x = coef, y = neg_log10_p,
                              color = Direction_pval, size = abs(coef))) +
    geom_point(alpha = 0.72) +
    geom_hline(yintercept = -log10(abundance_p_threshold),
               linetype = "dashed", color = "grey35", linewidth = 0.6) +
    geom_vline(xintercept = c(-abundance_coef_threshold,
                              abundance_coef_threshold),
               linetype = "dashed", color = "grey35", linewidth = 0.6) +
    scale_color_manual(values = abundance_volcano_colors) +
    scale_size_continuous(range = c(1.5, 4.5), guide = "none") +
    annotate("text",
             x = x_range[2] * 0.95, y = -log10(abundance_p_threshold),
             label    = sprintf("p = %.2f", abundance_p_threshold),
             hjust = 1, vjust = -0.5, size = 3.2,
             color = "grey35", fontface = "italic") +
    annotate("text",
             x = x_range[1], y = y_max_p,
             label    = sprintf("n = %d", n_down_p),
             hjust = 0, vjust = 1, size = 4,
             color = unname(group_colors[abundance_ref_group]),
             fontface = "bold") +
    annotate("text",
             x = x_range[2], y = y_max_p,
             label    = sprintf("n = %d", n_up_p),
             hjust = 1, vjust = 1, size = 4,
             color = unname(group_colors[abundance_comp_group]),
             fontface = "bold") +
    labs(
      title = sprintf("Volcano Plot \u2014 %s (p-value)", ds_name),
      subtitle = sprintf(
        "%s vs %s  \u00B7  MaAsLin2: TSS+LOG  \u00B7  uncorrected p < %.2f",
        abundance_comp_group, abundance_ref_group, abundance_p_threshold),
      x = sprintf("MaAsLin2 Coefficient\n\u2190 Higher in %s  |  Higher in %s \u2192",
                  abundance_ref_group, abundance_comp_group),
      y     = expression(-log[10](p-value)),
      color = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(color = "grey40", size = 9.5),
          legend.position = "bottom",
          legend.text     = element_text(size = 10))
  
  if (nrow(top_labels_p) > 0) {
    p_volc_p <- p_volc_p +
      geom_text_repel(data          = top_labels_p,
                      aes(label     = feature_clean),
                      size          = 2.7,
                      max.overlaps  = 30,
                      segment.color = "grey50",
                      segment.size  = 0.3,
                      box.padding   = 0.45,
                      show.legend   = FALSE)
  }
  
  # ═════════════════════════════════════════════════════════════════════════════
  #  EFFECT-SIZE BAR PLOT
  # ═════════════════════════════════════════════════════════════════════════════
  
  half_bar <- ceiling(abundance_n_top_bar / 2)
  
  top_pos <- res %>%
    filter(qval < abundance_q_threshold, coef > 0) %>%
    arrange(qval, desc(coef)) %>%
    head(half_bar)
  
  top_neg <- res %>%
    filter(qval < abundance_q_threshold, coef < 0) %>%
    arrange(qval, coef) %>%
    head(half_bar)
  
  top_bar_data <- bind_rows(top_neg, top_pos)
  
  p_bar   <- NULL
  p_panel <- NULL
  
  if (nrow(top_bar_data) > 0) {
    
    top_bar_data <- top_bar_data %>%
      mutate(feature_clean = forcats::fct_reorder(feature_clean, coef))
    
    p_bar <- ggplot(top_bar_data,
                    aes(x = coef, y = feature_clean, fill = Direction)) +
      geom_col(alpha = 0.85, width = 0.75) +
      geom_errorbarh(
        aes(xmin = coef - stderr, xmax = coef + stderr),
        height = 0.3, color = "grey25", linewidth = 0.5
      ) +
      geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
      scale_fill_manual(values = abundance_volcano_colors) +
      labs(
        title = sprintf("Top %d DA Features \u2014 %s",
                        nrow(top_bar_data), ds_name),
        subtitle = sprintf(
          "MaAsLin2  \u00B7  q < %.2f  \u00B7  error bars = \u00B11 SE",
          abundance_q_threshold),
        x = sprintf("Coefficient\n\u2190 %s  |  %s \u2192",
                    abundance_ref_group, abundance_comp_group),
        y = NULL, fill = NULL
      ) +
      theme_bw(base_size = 12) +
      theme(plot.title      = element_text(face = "bold"),
            plot.subtitle   = element_text(color = "grey40", size = 9.5),
            legend.position = "bottom",
            axis.text.y     = element_text(size = 8))
    
    # ── Painel combinado: volcano + barplot ──
    p_panel <- p_volc_q + p_bar +
      plot_layout(widths = c(1.3, 1)) +
      plot_annotation(
        title = sprintf("Differential Abundance \u2014 %s \u2014 MaAsLin2",
                        ds_name),
        theme = theme(plot.title = element_text(face = "bold", size = 15))
      )
    
  } else {
    cat("  AVISO: Nenhuma feature significativa (q < 0.25) para barplot.\n")
  }
  
  # ═════════════════════════════════════════════════════════════════════════════
  #  HEATMAP DE FEATURES DIFERENCIALMENTE ABUNDANTES (opcional)
  # ═════════════════════════════════════════════════════════════════════════════
  
  ht_da <- NULL
  
  if (n_sig_q25 >= 5) {
    
    cat(sprintf("\n  Gerando heatmap com %d features significativas...\n",
                n_sig_q25))
    
    sig_feats         <- res %>%
      filter(qval < abundance_q_threshold) %>% pull(feature)
    sig_feats_present <- intersect(sig_feats, colnames(features_df))
    
    if (length(sig_feats_present) >= 5) {
      
      feat_mat   <- as.matrix(features_df[, sig_feats_present, drop = FALSE])
      rel_ht     <- sweep(feat_mat, 1, rowSums(as.matrix(features_df)), "/")
      rel_ht[is.nan(rel_ht)] <- 0
      
      # Z-score por feature e transpor para features × samples
      hm_z <- t(scale(rel_ht))
      hm_z <- pmax(pmin(hm_z, 3), -3)
      hm_z[is.nan(hm_z)] <- 0
      
      da_heatmap_nfeats[[ds_name]] <- nrow(hm_z)
      
      hc_cols <- tryCatch(
        hclust(dist(t(hm_z)), method = "ward.D2"),
        error = function(e) FALSE)
      hc_rows <- tryCatch(
        hclust(dist(hm_z), method = "ward.D2"),
        error = function(e) FALSE)
      
      cf_ht <- colorRamp2(
        c(-3, -1.5, 0, 1.5, 3),
        c("#2166AC", "#92C5DE", "white", "#F4A582", "#D6604D"))
      
      sample_grps <- meta_df$Group[match(colnames(hm_z), rownames(meta_df))]
      ca_ht <- HeatmapAnnotation(
        Group = sample_grps,
        col   = list(Group = group_palette),
        annotation_name_gp = gpar(fontsize = 11, fontface = "bold"))
      
      rfs <- ifelse(nrow(hm_z) > 40, 5,
                    ifelse(nrow(hm_z) > 25, 6.5, 8))
      
      ht_da <- Heatmap(
        hm_z,
        name            = "Z-score",
        col             = cf_ht,
        top_annotation  = ca_ht,
        cluster_rows    = hc_rows,
        cluster_columns = hc_cols,
        show_column_names = TRUE,
        column_names_gp   = gpar(fontsize = 7),
        show_row_names    = TRUE,
        row_names_gp      = gpar(fontsize = rfs),
        row_names_side    = "left",
        row_dend_side     = "right",
        column_title = sprintf(
          "%s \u2014 %d Significant DA Features (q < %.2f)",
          ds_name, nrow(hm_z), abundance_q_threshold),
        column_title_gp = gpar(fontsize = 12, fontface = "bold"),
        heatmap_legend_param = list(
          title    = "Z-score",
          at       = c(-3, -1.5, 0, 1.5, 3),
          title_gp = gpar(fontface = "bold")),
        border = TRUE
      )
      draw(ht_da, merge_legend = TRUE)
    }
    
  } else if (n_sig_q25 > 0 && n_sig_q25 < 5) {
    cat(sprintf("  INFO: %d features significativas (min 5 para heatmap).\n",
                n_sig_q25))
  }
  
  # ── Armazenar ──────────────────────────────────────────────────────────────
  
  da_results_all[[ds_name]] <- res
  da_plots_all[[ds_name]]   <- list(volcano_q = p_volc_q,
                                    volcano_p = p_volc_p,
                                    bar       = p_bar,
                                    panel     = p_panel)
  da_heatmaps_all[[ds_name]] <- ht_da
  
  da_summary_all[[ds_name]] <- data.frame(
    Dataset             = ds_name,
    N_features_tested   = n_total,
    N_sig_q25           = n_sig_q25,
    N_sig_q10           = n_sig_q10,
    N_sig_q05           = n_sig_q05,
    N_sig_p05           = n_sig_p05,
    N_enriched_comp_q25 = n_up_q,
    N_enriched_ref_q25  = n_down_q,
    N_enriched_comp_p05 = n_up_p,
    N_enriched_ref_p05  = n_down_p,
    Reference_group     = abundance_ref_group,
    Comparison_group    = abundance_comp_group,
    stringsAsFactors    = FALSE
  )
  
  # Exibir gráficos
  print(p_volc_q)
  print(p_volc_p)
  if (!is.null(p_bar))   print(p_bar)
  if (!is.null(p_panel)) print(p_panel)
  
  cat(sprintf("\n  OK: %s completa!\n\n", ds_name))
}

# ═══════════════════════════════════════════════════════════════════════════════
# 11.6  RESUMO CONSOLIDADO
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n", rep("=", 70), "\n")
cat("=== RESUMO CONSOLIDADO \u2014 ABUNDANCIA DIFERENCIAL ===\n")
cat(rep("=", 70), "\n\n")

all_da_summary <- do.call(rbind, da_summary_all)

if (!is.null(all_da_summary) && nrow(all_da_summary) > 0) {
  print(all_da_summary)
} else {
  cat("  Nenhum resultado disponivel.\n")
}

# ── Gráfico panorâmico: nº de features DA por dataset ────────────────────────

p_overview <- NULL

if (!is.null(all_da_summary) && nrow(all_da_summary) > 0) {
  
  overview_long <- all_da_summary %>%
    select(Dataset, N_enriched_comp_q25, N_enriched_ref_q25) %>%
    rename(
      !!paste("Enriched in", abundance_comp_group) := N_enriched_comp_q25,
      !!paste("Enriched in", abundance_ref_group)  := N_enriched_ref_q25
    ) %>%
    tidyr::pivot_longer(-Dataset,
                        names_to  = "Direction",
                        values_to = "Count")
  
  p_overview <- ggplot(overview_long,
                       aes(x = reorder(Dataset, -Count),
                           y = Count, fill = Direction)) +
    geom_col(position = "dodge", alpha = 0.85, width = 0.7) +
    geom_text(aes(label = Count),
              position = position_dodge(width = 0.7),
              vjust = -0.3, size = 3.8) +
    scale_fill_manual(values = c(
      setNames(group_colors[abundance_comp_group],
               paste("Enriched in", abundance_comp_group)),
      setNames(group_colors[abundance_ref_group],
               paste("Enriched in", abundance_ref_group))
    )) +
    labs(
      title    = "Differential Abundance Overview \u2014 All Datasets",
      subtitle = sprintf(
        "MaAsLin2  \u00B7  q < %.2f  \u00B7  %s vs %s",
        abundance_q_threshold, abundance_comp_group, abundance_ref_group),
      x = "Dataset", y = "Number of DA Features", fill = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(plot.title      = element_text(face = "bold"),
          plot.subtitle   = element_text(color = "grey40", size = 10),
          legend.position = "bottom",
          axis.text.x     = element_text(angle = 25, hjust = 1))
  
  print(p_overview)
}

# ═══════════════════════════════════════════════════════════════════════════════
# 11.7  EXPORTAÇÃO — PLOTS
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n\nExportando resultados de abundancia diferencial...\n\n")

for (ds_name in names(da_results_all)) {
  
  ds_dir <- file.path(abundance_output_root, ds_name)
  dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)
  
  plots <- da_plots_all[[ds_name]]
  
  # Volcano q-value
  save_gg(plots$volcano_q,
          file.path(ds_dir, paste0(ds_name, "_volcano_qvalue")),
          width = 11, height = 8.5)
  
  # Volcano p-value
  save_gg(plots$volcano_p,
          file.path(ds_dir, paste0(ds_name, "_volcano_pvalue")),
          width = 11, height = 8.5)
  
  # Barplot (altura adaptável)
  if (!is.null(plots$bar)) {
    n_bar_feats <- nrow(
      da_results_all[[ds_name]] %>% filter(qval < abundance_q_threshold)
    )
    bar_h <- max(7, min(16, 4 + min(n_bar_feats, abundance_n_top_bar) * 0.35))
    save_gg(plots$bar,
            file.path(ds_dir, paste0(ds_name, "_barplot_effect_size")),
            width = 10, height = bar_h)
  }
  
  # Painel combinado
  if (!is.null(plots$panel)) {
    save_gg(plots$panel,
            file.path(ds_dir, paste0(ds_name, "_panel_volcano_bar")),
            width = 20, height = 10)
  }
  
  # Heatmap
  if (!is.null(da_heatmaps_all[[ds_name]])) {
    tryCatch({
      n_hf  <- ifelse(!is.null(da_heatmap_nfeats[[ds_name]]),
                      da_heatmap_nfeats[[ds_name]], 30)
      ht_h  <- max(8, min(18, 5 + n_hf * 0.25))
      png(file.path(ds_dir, paste0(ds_name, "_heatmap_DA_features.png")),
          width = 14, height = ht_h, units = "in", res = 300, bg = "white")
      draw(da_heatmaps_all[[ds_name]], merge_legend = TRUE)
      dev.off()
      cat(sprintf("    Salvo: %s_heatmap_DA_features.png\n", ds_name))
    }, error = function(e) {
      cat(sprintf("    AVISO: Heatmap export falhou para %s\n", ds_name))
    })
  }
}

# Overview plot
if (!is.null(p_overview)) {
  save_gg(p_overview,
          file.path(abundance_output_root, "overview_DA_all_datasets"),
          width = 12, height = 7)
}

# ═══════════════════════════════════════════════════════════════════════════════
# 11.8  EXPORTAÇÃO — TABELAS EXCEL
# ═══════════════════════════════════════════════════════════════════════════════

cat("\nExportando tabelas Excel...\n\n")

for (ds_name in names(da_results_all)) {
  
  ds_dir <- file.path(abundance_output_root, ds_name)
  
  wb_da <- createWorkbook()
  
  # All results
  export_cols <- intersect(
    c("Dataset", "feature", "feature_clean", "coef", "stderr",
      "pval", "qval", "Direction", "Direction_pval",
      "N", "N.not.0", "N.not.zero"),
    names(da_results_all[[ds_name]])
  )
  export_res <- da_results_all[[ds_name]] %>%
    select(all_of(export_cols)) %>%
    arrange(qval, pval)
  
  add_fmt_sheet(wb_da, "All_Results",
                as.data.frame(export_res), hdr_color = "#C0504D")
  
  # Significant (q-value)
  sig_q <- export_res %>% filter(qval < abundance_q_threshold)
  if (nrow(sig_q) > 0)
    add_fmt_sheet(wb_da, "Significant_qval",
                  as.data.frame(sig_q), hdr_color = "#9B2335")
  
  # Significant (p-value)
  sig_p <- export_res %>% filter(pval < abundance_p_threshold)
  if (nrow(sig_p) > 0)
    add_fmt_sheet(wb_da, "Significant_pval",
                  as.data.frame(sig_p), hdr_color = "#9B2335")
  
  # Summary
  add_fmt_sheet(wb_da, "Summary",
                da_summary_all[[ds_name]], hdr_color = "#7030A0")
  
  saveWorkbook(wb_da,
               file.path(ds_dir, paste0(ds_name, "_differential_abundance.xlsx")),
               overwrite = TRUE)
  cat(sprintf("    %s/%s_differential_abundance.xlsx\n", ds_name, ds_name))
}

# ── Master consolidado ────────────────────────────────────────────────────────

wm_abundance <- createWorkbook()

add_fmt_sheet(wm_abundance, "Summary",
              as.data.frame(all_da_summary), hdr_color = "#C0504D")

# All significant q-value
all_sig_q <- do.call(rbind, lapply(names(da_results_all), function(ds) {
  r <- da_results_all[[ds]] %>% filter(qval < abundance_q_threshold)
  if (nrow(r) == 0) return(NULL)
  r %>%
    select(any_of(c("Dataset", "feature", "feature_clean", "coef",
                    "stderr", "pval", "qval", "Direction"))) %>%
    arrange(qval)
}))
if (!is.null(all_sig_q) && nrow(all_sig_q) > 0)
  add_fmt_sheet(wm_abundance, "All_Significant_q",
                as.data.frame(all_sig_q), hdr_color = "#9B2335")

# All significant p-value
all_sig_p <- do.call(rbind, lapply(names(da_results_all), function(ds) {
  r <- da_results_all[[ds]] %>% filter(pval < abundance_p_threshold)
  if (nrow(r) == 0) return(NULL)
  r %>%
    select(any_of(c("Dataset", "feature", "feature_clean", "coef",
                    "stderr", "pval", "qval", "Direction_pval"))) %>%
    arrange(pval)
}))
if (!is.null(all_sig_p) && nrow(all_sig_p) > 0)
  add_fmt_sheet(wm_abundance, "All_Significant_p",
                as.data.frame(all_sig_p), hdr_color = "#9B2335")

# All results combined
all_res_combined <- do.call(rbind, lapply(names(da_results_all), function(ds) {
  da_results_all[[ds]] %>%
    select(any_of(c("Dataset", "feature", "feature_clean", "coef",
                    "stderr", "pval", "qval", "Direction"))) %>%
    arrange(qval)
}))
add_fmt_sheet(wm_abundance, "All_Results",
              as.data.frame(all_res_combined), hdr_color = "#C0504D")

saveWorkbook(wm_abundance,
             file.path(abundance_output_root,
                       "MASTER_abundance_all_datasets.xlsx"),
             overwrite = TRUE)
cat("\n    MASTER_abundance_all_datasets.xlsx\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 11.9  MENSAGEM FINAL
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n", rep("=", 70), "\n")
cat("  FASE 11 CONCLUIDA \u2014 ANALISE DE ABUNDANCIA DIFERENCIAL\n")
cat(rep("=", 70), "\n\n")
cat(sprintf("  Pasta: %s/\n\n", abundance_output_root))

for (ds_name in names(da_results_all)) {
  n_sig <- da_summary_all[[ds_name]]$N_sig_q25
  cat(sprintf("    %s/  (plots + tabela | %d sig q<0.25)\n", ds_name, n_sig))
}

cat("\n    overview_DA_all_datasets.png\n")
cat("    MASTER_abundance_all_datasets.xlsx\n\n")
cat(rep("=", 70), "\n")
cat("PIPELINE FASES 1-11 CONCLUIDO COM SUCESSO!\n")
cat(rep("=", 70), "\n")
# ═══════════════════════════════════════════════════════════════════════════════
# FASE 12 — EXPORT
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n", rep("=",70), "\n=== FASE 12: EXPORT ===\n", rep("=",70), "\n\n")

for (ds in names(datasets_grouped_list)) {
  for (d in c(file.path(output_root,ds,"plots","pca"),
              file.path(output_root,ds,"plots","alpha_diversity"),
              file.path(output_root,ds,"plots","beta_diversity"),
              file.path(output_root,ds,"plots","differential_abundance"),
              file.path(output_root,ds,"tables")))
    dir.create(d, recursive=TRUE, showWarnings=FALSE)
}
dir.create(file.path(output_root,"consolidated"), recursive=TRUE, showWarnings=FALSE)

afs <- function(wb, sn, data, hc="#4472C4") {
  addWorksheet(wb, sn)
  hs<-createStyle(fontColour="#FFFFFF",fgFill=hc,fontName="Calibri",fontSize=11,textDecoration="bold",
                  halign="center",valign="center",border="TopBottomLeftRight",borderColour="#FFFFFF",wrapText=TRUE)
  bs<-createStyle(fontName="Calibri",fontSize=10,border="TopBottomLeftRight",borderColour="#BFBFBF")
  as3<-createStyle(fontName="Calibri",fontSize=10,fgFill="#F2F7FF",border="TopBottomLeftRight",borderColour="#BFBFBF")
  writeData(wb,sn,data,headerStyle=hs)
  nr<-nrow(data); nc<-ncol(data)
  if(nr>0) { addStyle(wb,sn,bs,rows=2:(nr+1),cols=1:nc,gridExpand=TRUE)
    er<-seq(3,nr+1,by=2); if(length(er)>0) addStyle(wb,sn,as3,rows=er,cols=1:nc,gridExpand=TRUE) }
  setColWidths(wb,sn,1:nc,"auto"); freezePane(wb,sn,firstRow=TRUE); invisible(wb)
}

cat("Exportando plots...\n")
for (ds in names(datasets_grouped_list)) {
  cat(paste("  [",ds,"]\n"))
  # PCA
  dp<-file.path(output_root,ds,"plots","pca")
  if(!is.null(pca_results_grouped[[ds]])){
    sgg(pca_results_grouped[[ds]]$plot_pc12, file.path(dp,"pca_PC1_PC2"),9,7)
    if(!is.null(pca_results_grouped[[ds]]$plot_pc13)) sgg(pca_results_grouped[[ds]]$plot_pc13,file.path(dp,"pca_PC1_PC3"),9,7)
  }
  # Alpha
  da2<-file.path(output_root,ds,"plots","alpha_diversity")
  if(!is.null(alpha_plots_all[[ds]])){
    ap<-alpha_plots_all[[ds]]
    sgg(ap$depth,file.path(da2,"sequencing_depth"),13,6); sgg(ap$facet,file.path(da2,"alpha_all_metrics"),12,10)
    sgg(ap$panel,file.path(da2,"alpha_2x2_panel"),14,12); sgg(ap$observed,file.path(da2,"alpha_observed"),7,6)
    sgg(ap$chao1,file.path(da2,"alpha_chao1"),7,6); sgg(ap$shannon,file.path(da2,"alpha_shannon"),7,6)
    sgg(ap$simpson,file.path(da2,"alpha_simpson"),7,6)
  }
  # Beta
  db2<-file.path(output_root,ds,"plots","beta_diversity")
  if(!is.null(beta_plots_all[[ds]])){
    bp<-beta_plots_all[[ds]]
    for(dn in names(bp$pcoa_plots)) sgg(bp$pcoa_plots[[dn]],file.path(db2,paste0("pcoa_",gsub("[^a-zA-Z0-9]","_",dn))),9,7)
    sgg(bp$pcoa_panel,file.path(db2,"pcoa_panel"),16,7
