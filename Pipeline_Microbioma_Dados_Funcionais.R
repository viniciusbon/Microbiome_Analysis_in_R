###############################################################################
#                    CONFIGURAÇÃO DO USUÁRIO
###############################################################################

data_dir           <- "C:/Users/dti-/Desktop/Padronizados e prontos para analise/HAV2112"
sample_id_col_name <- "Sample_ID"
group_col_name     <- "group"
selected_groups    <- c("control", "symph900")
age_col_name       <- "age"
selected_ages      <- NULL
group_colors       <- c("#E31A1C", "#1F78B4")
names(group_colors) <- selected_groups

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
# FASE 1 — IMPORTAÇÃO
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n", rep("=",60), "\n=== FASE 1: IMPORTACAO ===\n", rep("=",60), "\n")
cat("Diretorio:", data_dir, "\n\n")

patterns <- c(AMR="AMR_counts", CAZy="CAZy_counts", COGs="COGs_count",
              EC="EC_count", KO="KO_counts", VF="VF_counts",
              TAX="tax_species_taxonomic_profiles", meta="meta")

data_storage <- list()
for (key in names(patterns)) {
  ff <- list.files(data_dir, pattern=patterns[key], full.names=TRUE, ignore.case=TRUE)
  if (length(ff)==0) { warning(paste("Not found:",key)); next }
  if (length(ff)>1) ff <- ff[1]
  data_storage[[key]] <- if (grepl("\\.xlsx$",ff,ignore.case=TRUE)) read_excel(ff) else
    read_delim(ff, delim="\t", show_col_types=FALSE)
  message(paste("Imported:", key, "-", basename(ff)))
}

df_ko<-data_storage$KO; df_AMR<-data_storage$AMR; df_VF<-data_storage$VF
df_CAZy<-data_storage$CAZy; df_EC<-data_storage$EC; df_COGs<-data_storage$COGs
df_TAX<-data_storage$TAX; df_meta<-data_storage$meta

for (col in c("sample","Sample","sample_id","SampleID","sample_ID","Sample_ID")) {
  if (col %in% colnames(df_meta) && col != sample_id_col_name) {
    colnames(df_meta)[colnames(df_meta)==col] <- sample_id_col_name
    cat(paste("Renomeada:", col, "->", sample_id_col_name, "\n")); break
  } else if (col %in% colnames(df_meta) && col==sample_id_col_name) {
    cat(paste(sample_id_col_name, "OK\n")); break
  }
}

cat("\n=== VERIFICACAO ===\n")
for (nm in c("KO","AMR","VF","CAZy","EC","COGs","TAX","META")) {
  obj <- switch(nm, KO=df_ko, AMR=df_AMR, VF=df_VF, CAZy=df_CAZy,
                EC=df_EC, COGs=df_COGs, TAX=df_TAX, META=df_meta)
  if (!is.null(obj)) cat(paste(nm,":",nrow(obj),"x",ncol(obj),"\n")) else cat(paste(nm,": NULL\n"))
}

# ═══════════════════════════════════════════════════════════════════════════════
# FASE 2 — TRANSPOSIÇÃO
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== FASE 2: TRANSPOSICAO ===\n")

transpose_fd <- function(data, label) {
  fids <- data[[1]]; sn <- names(data)[-1]; cm <- data[,-1]
  td <- as.data.frame(t(cm)); names(td) <- fids
  td[[sample_id_col_name]] <- sn
  td[, c(sample_id_col_name, setdiff(names(td), sample_id_col_name))]
}

for (nm in c("df_ko","df_CAZy","df_COGs","df_EC")) {
  obj <- get(nm, envir=.GlobalEnv)
  if (!is.null(obj)) {
    tr <- transpose_fd(obj, nm); assign(nm, tr, envir=.GlobalEnv)
    cat(paste("OK",nm,":",nrow(tr),"x",ncol(tr),"\n"))
  }
}
df_ko<-get("df_ko",.GlobalEnv); df_CAZy<-get("df_CAZy",.GlobalEnv)
df_COGs<-get("df_COGs",.GlobalEnv); df_EC<-get("df_EC",.GlobalEnv)

# FASE 2A — TAX
cat("\n=== FASE 2A: TAX ===\n")
tax_taxonomy_table <- NULL
if (!is.null(df_TAX)) {
  tcf <- intersect(names(df_TAX), c("k","p","c","o","f","g","s","Kingdom","Phylum","Class","Order","Family","Genus","Species"))
  if (length(tcf)>0) {
    tax_taxonomy_table <- df_TAX[, tcf, drop=FALSE]
    if ("g" %in% tcf && "s" %in% tcf) sl <- paste(df_TAX[["g"]],df_TAX[["s"]],sep="|")
    else if ("s" %in% tcf) sl <- df_TAX[["s"]] else sl <- df_TAX[[tail(tcf,1)]]
    sl <- make.unique(as.character(sl)); cc <- setdiff(names(df_TAX), tcf)
    td <- as.data.frame(t(as.matrix(df_TAX[,cc,drop=FALSE]))); colnames(td) <- sl
    td[[sample_id_col_name]] <- rownames(td)
    td <- td[, c(sample_id_col_name, setdiff(names(td), sample_id_col_name))]
    rownames(td) <- NULL; df_TAX <- td
    cat(paste("  OK TAX:",nrow(df_TAX),"x",ncol(df_TAX),"\n"))
  } else { df_TAX <- transpose_fd(df_TAX,"tax") }
} else cat("  TAX nao encontrado\n")

# FASE 2B — PIVOT
cat("\n=== FASE 2B: PIVOT ===\n")
pivot_lw <- function(dataset, dsn, id_col) {
  if (is.null(dataset)) return(NULL)
  if (id_col %in% names(dataset)) {
    nr <- nrow(dataset); nu <- length(unique(dataset[[id_col]]))
    if (nu==nr) { cat(paste(" ",dsn,": largo OK\n")); return(dataset) }
    cat(paste(" ",dsn,": LONGO (",nr,"->",nu,")\n"))
  } else return(dataset)
  pm <- c(id_col,"group","age","treatment","study","tissue","age_group","sample")
  nm2 <- setdiff(names(dataset), pm)
  fc <- NULL; for (cl in nm2) if (is.character(dataset[[cl]])||is.factor(dataset[[cl]])) { fc<-cl; break }
  if (is.null(fc)) fc <- nm2[1]
  nc2 <- nm2[sapply(dataset[,nm2,drop=FALSE], is.numeric)]
  ccands <- c("Depth","depth","Hits","hits","Count","count","Score","score")
  vc <- NULL; for (cd in ccands) if (cd %in% nc2) { vc<-cd; break }
  if (is.null(vc)) {
    pd <- dataset[,c(id_col,fc),drop=FALSE]; pd$cv <- 1
    ag <- aggregate(as.formula(paste("cv~",id_col,"+",fc)), data=pd, FUN=sum)
    wd <- reshape(ag, idvar=id_col, timevar=fc, direction="wide", v.names="cv")
    names(wd) <- gsub("^cv\\.","",names(wd))
  } else {
    pd <- dataset[,c(id_col,fc,vc),drop=FALSE]
    ag <- aggregate(as.formula(paste(vc,"~",id_col,"+",fc)), data=pd, FUN=sum)
    wd <- reshape(ag, idvar=id_col, timevar=fc, direction="wide", v.names=vc)
    pf <- paste0(vc,"."); names(wd) <- gsub(paste0("^",gsub("([.])","\\\\\\1",pf)),"",names(wd))
  }
  wd[is.na(wd)] <- 0; attr(wd,"reshapeWide") <- NULL; rownames(wd) <- NULL
  cat(paste("   Resultado:",nrow(wd),"x",ncol(wd),"\n")); wd
}

df_AMR<-pivot_lw(df_AMR,"AMR",sample_id_col_name); df_VF<-pivot_lw(df_VF,"VF",sample_id_col_name)
df_ko<-pivot_lw(df_ko,"KO",sample_id_col_name); df_CAZy<-pivot_lw(df_CAZy,"CAZy",sample_id_col_name)
df_EC<-pivot_lw(df_EC,"EC",sample_id_col_name); df_COGs<-pivot_lw(df_COGs,"COGs",sample_id_col_name)
df_TAX<-pivot_lw(df_TAX,"TAX",sample_id_col_name)

# ═══════════════════════════════════════════════════════════════════════════════
# FASE 3 — PADRONIZAÇÃO NOMES
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== FASE 3: PADRONIZACAO ===\n")
stdn <- function(x) x %>% str_replace_all("\\.","-") %>% str_replace("_$","")
if (!is.null(df_meta) && sample_id_col_name %in% names(df_meta))
  df_meta[[sample_id_col_name]] <- stdn(df_meta[[sample_id_col_name]])

# ═══════════════════════════════════════════════════════════════════════════════
# FASE 4 — FILTRAGEM
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== FASE 4: FILTRAGEM ===\n")
ri <- df_meta[[sample_id_col_name]]; cat("Ref:",length(ri),"\n")
flt <- function(ds,nm,r) {
  if (is.null(ds)||!sample_id_col_name %in% names(ds)) return(NULL)
  fd <- ds[ds[[sample_id_col_name]] %in% r,]
  cat(paste(" ",nm,":",nrow(ds),"->",nrow(fd),"\n")); fd
}
df_ko<-flt(df_ko,"KO",ri); df_AMR<-flt(df_AMR,"AMR",ri); df_VF<-flt(df_VF,"VF",ri)
df_CAZy<-flt(df_CAZy,"CAZy",ri); df_EC<-flt(df_EC,"EC",ri)
df_COGs<-flt(df_COGs,"COGs",ri); df_TAX<-flt(df_TAX,"TAX",ri)

# ═══════════════════════════════════════════════════════════════════════════════
# FASE 5 — JOIN METADATA
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== FASE 5: JOIN ===\n")
rmc <- c(sample_id_col_name, group_col_name, age_col_name)
miss <- setdiff(rmc, names(df_meta))
if (length(miss)>0) stop(paste("ERRO: Faltam:", paste(miss,collapse=", "),
                               "\nDisponiveis:", paste(names(df_meta),collapse=", ")))
msub <- df_meta[, rmc, drop=FALSE]

dl <- list(df_ko=df_ko, df_AMR=df_AMR, df_VF=df_VF, df_CAZy=df_CAZy,
           df_EC=df_EC, df_COGs=df_COGs, df_TAX=df_TAX)
dj <- lapply(names(dl), function(nm) {
  d<-dl[[nm]]; if(is.null(d)) return(NULL)
  m<-merge(d,msub,by=sample_id_col_name,all.x=TRUE)
  co<-c(sample_id_col_name,group_col_name,age_col_name,
        setdiff(names(m),c(sample_id_col_name,group_col_name,age_col_name)))
  m[,co]
})
names(dj)<-names(dl); list2env(dj,envir=.GlobalEnv)
df_ko<-dj$df_ko; df_AMR<-dj$df_AMR; df_VF<-dj$df_VF; df_CAZy<-dj$df_CAZy
df_EC<-dj$df_EC; df_COGs<-dj$df_COGs; df_TAX<-dj$df_TAX

# ═══════════════════════════════════════════════════════════════════════════════
# FASE 6 — FILTRO GRUPOS
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== FASE 6: FILTRO GRUPOS ===\n")
cat("Grupos:", paste(selected_groups,collapse=", "), "\n")
dff <- list(df_ko=df_ko, df_AMR=df_AMR, df_VF=df_VF, df_CAZy=df_CAZy,
            df_EC=df_EC, df_COGs=df_COGs, df_TAX=df_TAX)
dg <- lapply(names(dff), function(nm) {
  d<-dff[[nm]]; if(is.null(d)||!group_col_name %in% names(d)) return(NULL)
  fd<-d[d[[group_col_name]] %in% selected_groups,]
  if (!is.null(selected_ages)) fd<-fd[fd[[age_col_name]] %in% selected_ages,]
  gt<-table(fd[[group_col_name]])
  cat(paste(" ",nm,":",nrow(d),"->",nrow(fd),"|",paste(names(gt),gt,sep="=",collapse=", "),"\n")); fd
})
names(dg)<-paste0(names(dff),"_grouped"); list2env(dg,envir=.GlobalEnv)
df_ko_grouped<-dg$df_ko_grouped; df_AMR_grouped<-dg$df_AMR_grouped
df_VF_grouped<-dg$df_VF_grouped; df_CAZy_grouped<-dg$df_CAZy_grouped
df_EC_grouped<-dg$df_EC_grouped; df_COGs_grouped<-dg$df_COGs_grouped
df_TAX_grouped<-dg$df_TAX_grouped
# ═══════════════════════════════════════════════════════════════════════════════
# FUNÇÕES AUXILIARES GLOBAIS
# ═══════════════════════════════════════════════════════════════════════════════
meta_cols_internal <- c(sample_id_col_name, group_col_name, age_col_name)

ecm <- function(ds) {
  fc<-setdiff(names(ds),meta_cols_internal)
  m<-as.data.frame(lapply(ds[,fc],function(x) as.numeric(as.character(x))))
  m[is.na(m)]<-0; rownames(m)<-ds[[sample_id_col_name]]; as.matrix(m)
}
egv <- function(ds) { g<-ds[[group_col_name]]; names(g)<-ds[[sample_id_col_name]]; g }
ss <- function(p) ifelse(p<.001,"***",ifelse(p<.01,"**",ifelse(p<.05,"*","ns")))
save_gg <- function(po,fp,w=10,h=7,d=300) {
  ggsave(paste0(fp,".png"),plot=po,width=w,height=h,dpi=d,bg="white")
  cat(paste("    Salvo:",basename(fp),".png\n"))
}
detect_data_type <- function(mat) {
  vals<-mat[mat!=0]; if(length(vals)==0) return(list(type="empty",is_raw=FALSE))
  pct<-sum(vals==round(vals))/length(vals)*100
  list(type=ifelse(pct>95,"raw_counts","pre_normalized"), is_raw=pct>95)
}

# ═══════════════════════════════════════════════════════════════════════════════
# FASE 7 — PCA
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== FASE 7: PCA ===\n")

do_pca <- function(ds, dsn) {
  if (is.null(ds)||!group_col_name %in% names(ds)) return(NULL)
  cat(paste("PCA:",dsn,"\n"))
  fc<-setdiff(names(ds),meta_cols_internal); ms<-ds[,intersect(meta_cols_internal,names(ds)),drop=FALSE]
  fd<-as.data.frame(lapply(ds[,fc,drop=FALSE],function(x) as.numeric(as.character(x))))
  fd[is.na(fd)]<-0; vv<-apply(fd,2,var,na.rm=TRUE); vf<-!is.na(vv)&vv>0; fd<-fd[,vf,drop=FALSE]
  if(ncol(fd)<3) return(NULL)
  clr<-compositions::clr(fd+1e-6); pr<-prcomp(clr,center=TRUE,scale.=TRUE)
  np<-min(3,ncol(pr$x)); ve<-summary(pr)$importance["Proportion of Variance",1:np]*100
  pdf2<-as.data.frame(pr$x[,1:np,drop=FALSE])
  pdf2$Sample_ID<-ms[[sample_id_col_name]]; pdf2$Group<-ms[[group_col_name]]
  
  # Verificar se tem >= 4 amostras por grupo para elipse
  grp_n <- table(pdf2$Group)
  ellipse_data <- pdf2[pdf2$Group %in% names(grp_n[grp_n >= 4]), ]
  
  p12<-ggplot(pdf2,aes(PC1,PC2,color=Group,label=Sample_ID))+geom_point(size=3.5,alpha=.85)+
    geom_text_repel(size=2.5,max.overlaps=30,show.legend=FALSE)+
    scale_color_manual(values=group_colors)+
    labs(title=paste("PCA-",dsn,"PC1vPC2"),x=sprintf("PC1(%.1f%%)",ve[1]),y=sprintf("PC2(%.1f%%)",ve[2]))+
    theme_minimal(base_size=13)+theme(plot.title=element_text(face="bold"))
  if(nrow(ellipse_data)>=4) p12<-p12+stat_ellipse(data=ellipse_data,aes(group=Group),type="norm",linetype="dashed",alpha=.5)
  
  p13<-NULL
  if(np>=3) {
    p13<-ggplot(pdf2,aes(PC1,PC3,color=Group,label=Sample_ID))+geom_point(size=3.5,alpha=.85)+
      geom_text_repel(size=2.5,max.overlaps=30,show.legend=FALSE)+
      scale_color_manual(values=group_colors)+
      labs(title=paste("PCA-",dsn,"PC1vPC3"),x=sprintf("PC1(%.1f%%)",ve[1]),y=sprintf("PC3(%.1f%%)",ve[3]))+
      theme_minimal(base_size=13)+theme(plot.title=element_text(face="bold"))
    if(nrow(ellipse_data)>=4) p13<-p13+stat_ellipse(data=ellipse_data,aes(group=Group),type="norm",linetype="dashed",alpha=.5)
  }
  cat(paste("  OK:",paste(sprintf("%.1f%%",ve),collapse=", "),"\n"))
  list(pca_result=pr,pca_data=pdf2,variance_explained=ve,plot_pc12=p12,plot_pc13=p13,
       dataset_name=dsn,group_distribution=table(pdf2$Group))
}

dpca<-list(KO=df_ko_grouped,AMR=df_AMR_grouped,VF=df_VF_grouped,CAZy=df_CAZy_grouped,
           EC=df_EC_grouped,COGs=df_COGs_grouped,TAX=df_TAX_grouped)
pca_results_grouped<-list()
for(dn in names(dpca)) pca_results_grouped[[dn]]<-do_pca(dpca[[dn]],dn)
for(dn in names(pca_results_grouped)) if(!is.null(pca_results_grouped[[dn]])) print(pca_results_grouped[[dn]]$plot_pc12)

variance_summary<-do.call(rbind,lapply(names(pca_results_grouped),function(nm){
  r<-pca_results_grouped[[nm]]; if(is.null(r)) return(NULL)
  ve<-r$variance_explained; np<-length(ve)
  data.frame(Dataset=nm,PC1=round(ve[1],1),PC2=if(np>=2)round(ve[2],1)else NA,
             PC3=if(np>=3)round(ve[3],1)else NA,Total=round(sum(ve),1),stringsAsFactors=FALSE)
}))
for(g in selected_groups) variance_summary[[paste0(g,"_n")]]<-sapply(variance_summary$Dataset,function(d){
  gd<-pca_results_grouped[[d]]$group_distribution; if(g %in% names(gd)) as.integer(gd[g]) else 0L})
print(variance_summary)

# ═══════════════════════════════════════════════════════════════════════════════
# FASE 8 — ALPHA DIVERSITY
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n",rep("=",70),"\n=== FASE 8: ALPHA DIVERSITY ===\n",rep("=",70),"\n\n")

datasets_grouped_list<-list(KO=df_ko_grouped,AMR=df_AMR_grouped,VF=df_VF_grouped,
                            CAZy=df_CAZy_grouped,EC=df_EC_grouped,COGs=df_COGs_grouped,TAX=df_TAX_grouped)
datasets_grouped_list<-datasets_grouped_list[!sapply(datasets_grouped_list,is.null)]

output_root<-file.path(data_dir,"functional_diversity_results")
dir.create(output_root,recursive=TRUE,showWarnings=FALSE)

ng<-length(selected_groups)
if(ng<=2) group_palette<-setNames(c("#E31A1C","#1F78B4")[1:ng],selected_groups) else
  group_palette<-setNames(brewer.pal(min(ng,9),"Set1")[1:ng],selected_groups)
if(ng==2) comparisons_list<-list(selected_groups) else comparisons_list<-combn(selected_groups,2,simplify=FALSE)

alpha_results_all<-list(); alpha_stats_all<-list(); alpha_plots_all<-list(); rarefaction_info_all<-list()
data_type_info<-list()

for(ds in names(datasets_grouped_list)) {
  cat(rep("-",60),"\n"); cat(paste("ALPHA:",ds,"\n")); cat(rep("-",60),"\n")
  dat<-datasets_grouped_list[[ds]]; cm<-ecm(dat); gv<-egv(dat)
  dtype<-detect_data_type(cm); data_type_info[[ds]]<-dtype
  cat(sprintf("  Tipo: %s\n",dtype$type))
  cm[cm<0]<-0; cm<-cm[,colSums(cm)>0,drop=FALSE]
  sd2<-rowSums(cm); rd<-min(sd2)
  cat(paste("  Matriz:",nrow(cm),"x",ncol(cm),"| Depth:",round(rd),"\n"))
  
  rarefaction_info_all[[ds]]<-data.frame(Dataset=ds,Type=dtype$type,Min=round(rd),
                                         Max=round(max(sd2)),Med=round(median(sd2)),Samples=nrow(cm),Features=ncol(cm),stringsAsFactors=FALSE)
  
  ddf<-data.frame(Sample=names(sd2),Depth=as.numeric(sd2),Group=gv[names(sd2)],stringsAsFactors=FALSE)
  rare_label<-ifelse(dtype$is_raw,sprintf("Rarefied to %d",round(rd)),"Pre-normalized")
  pd<-ggplot(ddf,aes(reorder(Sample,Depth),Depth,fill=Group))+geom_col(alpha=.85)+
    scale_fill_manual(values=group_palette)+labs(title=paste(ds,"Depth"),x="Sample",y="Counts")+
    theme_minimal()+theme(axis.text.x=element_text(angle=90,hjust=1,size=7))
  if(dtype$is_raw) pd<-pd+geom_hline(yintercept=rd,color="darkred",linetype="dashed")
  
  if(dtype$is_raw) {
    cm_int<-round(cm); if(rd>=10){set.seed(42);cr<-vegan::rrarefy(cm_int,rd);cr<-cr[,colSums(cr)>0,drop=FALSE]
    } else cr<-cm_int
  } else cr<-cm
  
  or<-rowSums(cr>0)
  ch<-tryCatch({if(dtype$is_raw) t(vegan::estimateR(round(cr)))[,"S.chao1"] else or},error=function(e) or)
  sh<-vegan::diversity(cr,"shannon"); si<-vegan::diversity(cr,"simpson")
  adf<-data.frame(Sample=rownames(cr),Group=gv[rownames(cr)],Observed=as.numeric(or),
                  Chao1=as.numeric(ch),Shannon=as.numeric(sh),Simpson=as.numeric(si),stringsAsFactors=FALSE)
  alpha_results_all[[ds]]<-adf
  
  al<-tidyr::pivot_longer(adf,c("Observed","Chao1","Shannon","Simpson"),names_to="Metric",values_to="Value")
  al$Metric<-factor(al$Metric,levels=c("Observed","Chao1","Shannon","Simpson"))
  
  paf<-ggplot(al,aes(Group,Value,fill=Group,color=Group))+
    geom_violin(alpha=.35,trim=FALSE)+geom_boxplot(width=.18,alpha=.8,outlier.shape=NA,color="black")+
    geom_jitter(width=.08,size=1.8,alpha=.7,show.legend=FALSE)+
    stat_compare_means(method="wilcox.test",comparisons=comparisons_list,label="p.format",tip.length=.01,size=3.5)+
    scale_fill_manual(values=group_palette)+scale_color_manual(values=group_palette)+
    facet_wrap(~Metric,scales="free_y",nrow=2)+
    labs(title=paste("Alpha -",ds),subtitle=paste(rare_label,"| Wilcoxon"),x=NULL,y="Value")+
    theme_bw(base_size=13)+theme(plot.title=element_text(face="bold"),strip.text=element_text(face="bold"),
                                 legend.position="bottom",axis.text.x=element_text(angle=20,hjust=1))
  
  pam<-function(m,yl){d2<-al[al$Metric==m,]; ggplot(d2,aes(Group,Value,fill=Group,color=Group))+
    geom_violin(alpha=.3,trim=FALSE)+geom_boxplot(width=.2,alpha=.8,outlier.shape=NA,color="black")+
    geom_jitter(width=.08,size=2.2,alpha=.75,show.legend=FALSE)+
    stat_compare_means(method="wilcox.test",comparisons=comparisons_list,label="p.format",size=3.8)+
    scale_fill_manual(values=group_palette)+scale_color_manual(values=group_palette)+
    labs(title=m,x=NULL,y=yl)+theme_bw(base_size=13)+
    theme(plot.title=element_text(face="bold",hjust=.5),legend.position="none")}
  po<-pam("Observed","Richness"); pc<-pam("Chao1","Chao1"); ps<-pam("Shannon","Shannon"); pi2<-pam("Simpson","Simpson")
  pap<-(po|pc)/(ps|pi2)+plot_annotation(title=paste("Alpha -",ds),subtitle=rare_label)
  
  ml<-c("Observed","Chao1","Shannon","Simpson")
  kw<-do.call(rbind,lapply(ml,function(m){t2<-kruskal.test(adf[[m]]~factor(adf$Group))
  data.frame(Dataset=ds,Metric=m,Statistic=round(t2$statistic,4),p_value=round(t2$p.value,4),
             Sig=ss(t2$p.value),stringsAsFactors=FALSE)}))
  wx<-do.call(rbind,lapply(ml,function(m){do.call(rbind,lapply(comparisons_list,function(pr){
    va<-adf[adf$Group==pr[1],m]; vb<-adf[adf$Group==pr[2],m]; t2<-wilcox.test(va,vb,exact=FALSE)
    data.frame(Dataset=ds,Metric=m,Comparison=paste(pr,collapse=" vs "),W=round(t2$statistic,4),
               p_value=round(t2$p.value,4),Sig=ss(t2$p.value),stringsAsFactors=FALSE)}))}))
  de<-do.call(rbind,lapply(ml,function(m){do.call(rbind,lapply(selected_groups,function(g){
    v<-adf[adf$Group==g,m]
    data.frame(Dataset=ds,Group=g,Metric=m,N=length(v),Mean=round(mean(v),4),Median=round(median(v),4),
               SD=round(sd(v),4),Min=round(min(v),4),Max=round(max(v),4),stringsAsFactors=FALSE)}))}))
  
  alpha_stats_all[[ds]]<-list(kruskal_wallis=kw,wilcoxon=wx,descriptive=de)
  alpha_plots_all[[ds]]<-list(depth=pd,facet=paf,panel=pap,observed=po,chao1=pc,shannon=ps,simpson=pi2)
  cat("\n  KW:\n"); print(kw[,c("Metric","Statistic","p_value","Sig")])
  print(pd); print(paf); print(pap); cat(paste("\n  OK:",ds,"\n\n"))
}
all_kw<-do.call(rbind,lapply(alpha_stats_all,function(x) x$kruskal_wallis))
all_wx<-do.call(rbind,lapply(alpha_stats_all,function(x) x$wilcoxon))
all_desc<-do.call(rbind,lapply(alpha_stats_all,function(x) x$descriptive))

# ═══════════════════════════════════════════════════════════════════════════════
# FASE 9 — BETA DIVERSITY
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n",rep("=",70),"\n=== FASE 9: BETA DIVERSITY ===\n",rep("=",70),"\n\n")

beta_stats_all<-list(); beta_plots_all<-list()

for(ds in names(datasets_grouped_list)) {
  cat(rep("-",60),"\n"); cat(paste("BETA:",ds,"\n")); cat(rep("-",60),"\n")
  dat<-datasets_grouped_list[[ds]]; dtype<-data_type_info[[ds]]
  cm<-ecm(dat); gv<-egv(dat); cm[cm<0]<-0; cm<-cm[,colSums(cm)>0,drop=FALSE]
  rd<-min(rowSums(cm))
  
  if(dtype$is_raw){cm_int<-round(cm);if(rd>=10){set.seed(42);cr<-vegan::rrarefy(cm_int,rd);cr<-cr[,colSums(cr)>0,drop=FALSE]}else cr<-cm_int
  } else cr<-cm
  
  rm2<-sweep(cr,1,rowSums(cr),"/"); rm2[is.nan(rm2)]<-0
  mb<-data.frame(Sample=rownames(cr),Group=gv[rownames(cr)],stringsAsFactors=FALSE)
  beta_label<-ifelse(dtype$is_raw,sprintf("Rarefied %d",round(rd)),"Pre-normalized")
  cat(paste("  Matriz:",nrow(cr),"x",ncol(cr),"|",beta_label,"\n"))
  
  db<-vegdist(rm2,"bray"); dj2<-vegdist(if(dtype$is_raw) round(cr) else cr,"jaccard",binary=TRUE)
  dl2<-list("Bray-Curtis"=db,"Jaccard"=dj2)
  
  sbd<-list()
  for(dn in names(dl2)){
    set.seed(42);pm<-adonis2(dl2[[dn]]~Group,data=mb,permutations=999,by="margin")
    pR2<-round(pm$R2[1],3);pF2<-round(pm$F[1],3);pp<-pm$`Pr(>F)`[1]
    set.seed(42);an<-anosim(dl2[[dn]],mb$Group,permutations=999);aR<-round(an$statistic,3);ap<-an$signif
    at2<-sprintf("PERMANOVA: R2=%.3f p=%.4f %s\nANOSIM: R=%.3f p=%.4f %s",pR2,pp,ss(pp),aR,ap,ss(ap))
    cat(sprintf("  [%s] R2=%.3f p=%.4f%s | R=%.3f p=%.4f%s\n",dn,pR2,pp,ss(pp),aR,ap,ss(ap)))
    sbd[[dn]]<-list(permanova_R2=pR2,permanova_F=pF2,permanova_p=pp,permanova_sig=ss(pp),
                    anosim_R=aR,anosim_p=ap,anosim_sig=ss(ap),annot_text=at2)
  }
  
  bdr<-do.call(rbind,lapply(names(dl2),function(dn){bd<-betadisper(dl2[[dn]],mb$Group);set.seed(42)
  pd2<-permutest(bd,permutations=999)
  data.frame(Dataset=ds,Distance=dn,F_stat=round(pd2$tab[1,"F"],4),
             p_value=round(pd2$tab[1,"Pr(>F)"],4),stringsAsFactors=FALSE)}))
  
  # Verificar elipses
  grp_n_beta <- table(mb$Group)
  ellipse_ok <- all(grp_n_beta >= 4)
  
  pp_list<-list()
  for(dn in names(dl2)){
    pr2<-pcoa(dl2[[dn]]);ei<-pr2$values$Eigenvalues;vp<-round(ei[ei>0]/sum(ei[ei>0])*100,1)
    pd3<-as.data.frame(pr2$vectors[,1:2]);colnames(pd3)<-c("A1","A2")
    pd3$Sample<-rownames(pd3);pd3$Group<-mb$Group[match(pd3$Sample,mb$Sample)]
    p<-ggplot(pd3,aes(A1,A2,color=Group,label=Sample))+geom_point(size=3.5,alpha=.85)+
      geom_text_repel(size=2.5,max.overlaps=25,show.legend=FALSE)+
      scale_color_manual(values=group_palette)+
      annotate("label",x=Inf,y=-Inf,label=sbd[[dn]]$annot_text,hjust=1.03,vjust=-.15,
               size=2.8,fill="white",alpha=.88,family="mono")+
      labs(title=sprintf("PCoA-%s-%s",ds,dn),subtitle=beta_label,
           x=sprintf("Axis1(%.1f%%)",vp[1]),y=sprintf("Axis2(%.1f%%)",vp[2]))+
      theme_bw(base_size=13)+theme(plot.title=element_text(face="bold"))
    if(ellipse_ok) p<-p+stat_ellipse(aes(group=Group),type="norm",linetype="dashed",alpha=.6)
    pp_list[[dn]]<-p; print(p)
  }
  pcp<-wrap_plots(pp_list,ncol=2)+plot_layout(guides="collect")+
    plot_annotation(title=paste("PCoA-",ds))&theme(legend.position="bottom"); print(pcp)
  
  np_list<-list(); nst<-list()
  for(dn in names(dl2)){set.seed(42)
    nr2<-tryCatch(metaMDS(dl2[[dn]],k=2,trymax=100,trace=FALSE,autotransform=FALSE),error=function(e) NULL)
    if(is.null(nr2)) next; sv<-nr2$stress; nst[[dn]]<-sv
    sl2<-ifelse(sv<.1,"excelente",ifelse(sv<.2,"bom",ifelse(sv<.3,"aceitavel","ruim")))
    ns2<-as.data.frame(scores(nr2,display="sites")); ns2$Sample<-rownames(ns2)
    ns2$Group<-mb$Group[match(ns2$Sample,mb$Sample)]
    p<-ggplot(ns2,aes(NMDS1,NMDS2,color=Group,label=Sample))+geom_point(size=3.5,alpha=.85)+
      geom_text_repel(size=2.5,max.overlaps=25,show.legend=FALSE)+
      scale_color_manual(values=group_palette)+
      annotate("text",x=-Inf,y=Inf,label=sprintf("Stress:%.4f(%s)",sv,sl2),hjust=-.1,vjust=1.4,size=3.5,color="grey30",fontface="italic")+
      labs(title=sprintf("NMDS-%s-%s",ds,dn),x="NMDS1",y="NMDS2")+
      theme_bw(base_size=13)+theme(plot.title=element_text(face="bold"))
    if(ellipse_ok) p<-p+stat_ellipse(aes(group=Group),type="norm",linetype="dashed",alpha=.6)
    np_list[[dn]]<-p; print(p)
  }
  pnp<-NULL; if(length(np_list)>0){pnp<-wrap_plots(np_list,ncol=2)+plot_layout(guides="collect")+
    plot_annotation(title=paste("NMDS-",ds))&theme(legend.position="bottom"); print(pnp)}
  
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
  at3<-do.call(rbind,lapply(names(sbd),function(dn){s<-sbd[[dn]]
  data.frame(Dataset=ds,Distance=dn,R_stat=s$anosim_R,p_value=s$anosim_p,Sig=s$anosim_sig,stringsAsFactors=FALSE)}))
  nt2<-if(length(nst)>0) data.frame(Dataset=ds,Distance=names(nst),Stress=round(unlist(nst),4),stringsAsFactors=FALSE) else
    data.frame(Dataset=character(),Distance=character(),Stress=numeric(),stringsAsFactors=FALSE)
  
  beta_stats_all[[ds]]<-list(permanova=pt2,anosim=at3,betadisper=bdr,nmds_stress=nt2)
  beta_plots_all[[ds]]<-list(pcoa_plots=pp_list,pcoa_panel=pcp,nmds_plots=np_list,nmds_panel=pnp,heatmap=ht2)
  cat(paste("\n  OK:",ds,"\n\n"))
}
all_permanova<-do.call(rbind,lapply(beta_stats_all,function(x) x$permanova))
all_anosim<-do.call(rbind,lapply(beta_stats_all,function(x) x$anosim))
all_betadisp<-do.call(rbind,lapply(beta_stats_all,function(x) x$betadisper))
all_nmds_str<-do.call(rbind,lapply(beta_stats_all,function(x) x$nmds_stress))

# ═══════════════════════════════════════════════════════════════════════════════
# FASE 11 — MAASLIN2
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n",rep("=",70),"\n=== FASE 11: MAASLIN2 ===\n",rep("=",70),"\n\n")

q_threshold<-0.25; coef_threshold<-0.5; n_top_bar<-20
ref_group<-selected_groups[1]; comparison_group<-selected_groups[2]
cat(sprintf("Ref: %s | Comp: %s | q: %.2f\n\n",ref_group,comparison_group,q_threshold))

volcano_colors<-c(
  setNames(group_palette[comparison_group],paste("Enriched in",comparison_group)),
  setNames(group_palette[ref_group],paste("Enriched in",ref_group)),
  "Not significant"="grey72")

da_results_all<-list(); da_plots_all<-list(); da_stats_all<-list()

for(ds in names(datasets_grouped_list)){
  cat(rep("-",60),"\n"); cat(paste("DA:",ds,"\n")); cat(rep("-",60),"\n")
  dat<-datasets_grouped_list[[ds]]; dtype<-data_type_info[[ds]]
  fc<-setdiff(names(dat),meta_cols_internal)
  mf<-as.data.frame(lapply(dat[,fc],function(x) as.numeric(as.character(x))))
  mf[is.na(mf)]<-0; rownames(mf)<-dat[[sample_id_col_name]]; mf<-mf[,colSums(mf)>0,drop=FALSE]
  mm<-data.frame(Group=dat[[group_col_name]],stringsAsFactors=FALSE)
  rownames(mm)<-dat[[sample_id_col_name]]; mm$Group<-relevel(factor(mm$Group),ref=ref_group)
  
  if(dtype$is_raw){mn<-"TSS";mt<-"LOG";nl<-"TSS+LOG (raw)"} else {mn<-"NONE";mt<-"LOG";nl<-"NONE+LOG (pre-norm)"}
  cat(sprintf("  Input: %d x %d | %s\n",nrow(mf),ncol(mf),nl))
  
  mod<-file.path(output_root,ds,"maaslin2_output"); dir.create(mod,recursive=TRUE,showWarnings=FALSE)
  set.seed(42)
  fit<-tryCatch(Maaslin2(input_data=mf,input_metadata=mm,output=mod,fixed_effects="Group",
                         normalization=mn,transform=mt,analysis_method="LM",min_abundance=0,min_prevalence=0,
                         max_significance=0.25,correction="BH",standardize=FALSE,plot_heatmap=FALSE,
                         plot_scatter=FALSE,cores=1),error=function(e){cat(paste("  ERRO:",e$message,"\n"));NULL})
  if(is.null(fit)){cat("  Pulando.\n\n");next}
  
  mr<-fit$results; mr<-mr[mr$metadata=="Group",]; mr<-mr[order(mr$qval,mr$pval),]
  mr$feature_clean<-sub("^[^|]*\\|","",mr$feature)
  mr$qval_safe<-pmax(mr$qval,1e-10); mr$neg_log10_q<- -log10(mr$qval_safe)
  mr$Direction<-ifelse(mr$qval<q_threshold & mr$coef>0,paste("Enriched in",comparison_group),
                       ifelse(mr$qval<q_threshold & mr$coef<0,paste("Enriched in",ref_group),"Not significant"))
  
  n25<-sum(mr$qval<.25,na.rm=TRUE); n05<-sum(mr$qval<.05,na.rm=TRUE)
  nec<-sum(mr$qval<q_threshold & mr$coef>0,na.rm=TRUE)
  ner<-sum(mr$qval<q_threshold & mr$coef<0,na.rm=TRUE)
  cat(sprintf("  Tested:%d | q<.25:%d | Enr %s:%d | Enr %s:%d\n",
              nrow(mr),n25,comparison_group,nec,ref_group,ner))
  
  sf<-mr[mr$qval<q_threshold,]
  tl<-if(nrow(sf)>0) do.call(rbind,lapply(unique(sf$Direction),function(d){
    s2<-sf[sf$Direction==d,]; s2<-s2[order(s2$qval,-abs(s2$coef)),]; head(s2,12)})) else NULL
  
  pv<-ggplot(mr,aes(coef,neg_log10_q,color=Direction,size=abs(coef)))+
    geom_point(alpha=.72)+
    geom_hline(yintercept=-log10(q_threshold),linetype="dashed",color="grey35")+
    geom_vline(xintercept=c(-coef_threshold,coef_threshold),linetype="dashed",color="grey35")+
    scale_color_manual(values=volcano_colors)+scale_size_continuous(range=c(1.5,4.5),guide="none")+
    labs(title=paste("DA -",ds),subtitle=sprintf("%s vs %s | %s | q=%.2f",comparison_group,ref_group,nl,q_threshold),
         x=sprintf("Coefficient\n<- %s | %s ->",ref_group,comparison_group),
         y=expression(-log[10](q)),color=NULL)+
    theme_bw(base_size=13)+theme(plot.title=element_text(face="bold"),legend.position="bottom")
  if(!is.null(tl)&&nrow(tl)>0) pv<-pv+geom_text_repel(data=tl,aes(label=feature_clean),size=2.9,
                                                      max.overlaps=30,segment.color="grey50",box.padding=.45,show.legend=FALSE)
  if(nrow(mr)>0) pv<-pv+
    annotate("text",x=min(mr$coef,na.rm=TRUE),y=max(mr$neg_log10_q,na.rm=TRUE),
             label=sprintf("n=%d",ner),hjust=0,vjust=1,size=4,color=group_palette[ref_group],fontface="bold")+
    annotate("text",x=max(mr$coef,na.rm=TRUE),y=max(mr$neg_log10_q,na.rm=TRUE),
             label=sprintf("n=%d",nec),hjust=1,vjust=1,size=4,color=group_palette[comparison_group],fontface="bold")
  print(pv)
  
  pb<-NULL; ppnl<-NULL
  tp<-head(mr[mr$qval<q_threshold & mr$coef>0,],n_top_bar/2)
  tn<-head(mr[mr$qval<q_threshold & mr$coef<0,],n_top_bar/2)
  tb<-rbind(tn,tp)
  if(nrow(tb)>0){
    tb$feature_clean<-factor(tb$feature_clean,levels=tb$feature_clean[order(tb$coef)])
    pb<-ggplot(tb,aes(coef,feature_clean,fill=Direction))+geom_col(alpha=.85,width=.75)+
      geom_errorbar(aes(xmin=coef-stderr,xmax=coef+stderr),width=.3,color="grey25",orientation="y")+
      geom_vline(xintercept=0)+scale_fill_manual(values=volcano_colors)+
      labs(title=sprintf("Top %d DA - %s",nrow(tb),ds),subtitle=sprintf("%s | q<%.2f",nl,q_threshold),
           x=sprintf("Coef\n<- %s | %s ->",ref_group,comparison_group),y=NULL,fill=NULL)+
      theme_bw(base_size=12)+theme(plot.title=element_text(face="bold"),legend.position="bottom")
    print(pb)
    ppnl<-pv+pb+plot_layout(widths=c(1.4,1))+plot_annotation(title=paste("DA -",ds)); print(ppnl)
  }
  
  dex<-data.frame(Dataset=ds,Feature=mr$feature_clean,Coefficient=round(mr$coef,4),
                  SE=round(mr$stderr,4),p_value=round(mr$pval,6),q_value=round(mr$qval,6),
                  Direction=mr$Direction,Normalization=nl,stringsAsFactors=FALSE)
  if("N" %in% names(mr)) dex$N<-mr$N
  if("N.not.0" %in% names(mr)) dex$N_not_zero<-mr$N.not.0
  
  da_results_all[[ds]]<-dex
  da_stats_all[[ds]]<-data.frame(Dataset=ds,Tested=nrow(mr),Sig_q025=n25,Sig_q005=n05,
                                 Enr_comp=nec,Enr_ref=ner,Normalization=nl,stringsAsFactors=FALSE)
  da_plots_all[[ds]]<-list(volcano=pv,bar=pb,panel=ppnl)
  cat(paste("\n  OK:",ds,"\n\n"))
}

all_da_stats<-do.call(rbind,da_stats_all)
all_da_results<-do.call(rbind,da_results_all)
all_da_sig<-if(!is.null(all_da_results)&&nrow(all_da_results)>0) all_da_results[all_da_results$q_value<q_threshold,] else data.frame()
# ═══════════════════════════════════════════════════════════════════════════════
# FASE 12 — EXPORT UNIFICADO COMPLETO (AUDITADO)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n",rep("=",70),"\n=== FASE 12: EXPORT UNIFICADO ===\n",rep("=",70),"\n\n")

# 12.0 Criar TODAS as pastas
for(ds in names(datasets_grouped_list)){
  for(sub in c("plots/pca","plots/alpha_diversity","plots/beta_diversity",
               "plots/differential_abundance","tables"))
    dir.create(file.path(output_root,ds,sub),recursive=TRUE,showWarnings=FALSE)
}
dir.create(file.path(output_root,"consolidated"),recursive=TRUE,showWarnings=FALSE)
cat("Pastas criadas para:",paste(names(datasets_grouped_list),collapse=", "),"\n\n")

ec <- list(plots=0, tables=0, errors=0)  # export counter

add_fmt_sheet <- function(wb, sn, data, hdr_color="#4472C4") {
  addWorksheet(wb, sn)
  hs<-createStyle(fontColour="#FFFFFF",fgFill=hdr_color,fontName="Calibri",fontSize=11,
                  textDecoration="bold",halign="center",valign="center",border="TopBottomLeftRight",
                  borderColour="#FFFFFF",wrapText=TRUE)
  bs<-createStyle(fontName="Calibri",fontSize=10,border="TopBottomLeftRight",borderColour="#BFBFBF")
  as2<-createStyle(fontName="Calibri",fontSize=10,fgFill="#F2F7FF",border="TopBottomLeftRight",borderColour="#BFBFBF")
  writeData(wb,sn,data,headerStyle=hs)
  nr<-nrow(data); nc<-ncol(data)
  if(nr>0){addStyle(wb,sn,bs,rows=2:(nr+1),cols=1:nc,gridExpand=TRUE)
    er<-seq(3,nr+1,by=2); if(length(er)>0) addStyle(wb,sn,as2,rows=er,cols=1:nc,gridExpand=TRUE)}
  setColWidths(wb,sn,1:nc,"auto"); freezePane(wb,sn,firstRow=TRUE); invisible(wb)
}

# ─────────────────────────────────────────────────────────────────────────────
# 12.1 EXPORT POR DATASET
# ─────────────────────────────────────────────────────────────────────────────

for(ds in names(datasets_grouped_list)){
  cat(rep("-",50),"\n"); cat(sprintf("  EXPORTANDO: %s\n",ds)); cat(rep("-",50),"\n")
  
  # ── PCA (Fase 7): plot_pc12, plot_pc13 ────────────────────────────────────
  dp<-file.path(output_root,ds,"plots","pca")
  if(!is.null(pca_results_grouped[[ds]])){
    pr<-pca_results_grouped[[ds]]
    if(!is.null(pr$plot_pc12)){tryCatch({save_gg(pr$plot_pc12,file.path(dp,"pca_PC1_vs_PC2"),9,7);ec$plots<-ec$plots+1},error=function(e){cat(sprintf("    ERRO: %s\n",e$message));ec$errors<<-ec$errors+1})}
    if(!is.null(pr$plot_pc13)){tryCatch({save_gg(pr$plot_pc13,file.path(dp,"pca_PC1_vs_PC3"),9,7);ec$plots<-ec$plots+1},error=function(e){cat(sprintf("    ERRO: %s\n",e$message));ec$errors<<-ec$errors+1})}
  }
  
  # ── ALPHA (Fase 8): depth, facet, panel, observed, chao1, shannon, simpson ─
  da<-file.path(output_root,ds,"plots","alpha_diversity")
  if(!is.null(alpha_plots_all[[ds]])){
    ap<-alpha_plots_all[[ds]]
    for(item in list(
      list(p=ap$depth,   n="sequencing_depth",  w=13,h=6),
      list(p=ap$facet,   n="alpha_all_metrics", w=12,h=10),
      list(p=ap$panel,   n="alpha_2x2_panel",   w=14,h=12),
      list(p=ap$observed,n="alpha_observed",     w=7, h=6),
      list(p=ap$chao1,   n="alpha_chao1",        w=7, h=6),
      list(p=ap$shannon, n="alpha_shannon",      w=7, h=6),
      list(p=ap$simpson, n="alpha_simpson",      w=7, h=6)
    )){if(!is.null(item$p)){tryCatch({save_gg(item$p,file.path(da,item$n),item$w,item$h);ec$plots<-ec$plots+1},
                                     error=function(e){cat(sprintf("    ERRO %s: %s\n",item$n,e$message));ec$errors<<-ec$errors+1})}}
  }
  
  # ── BETA (Fase 9): pcoa×2, pcoa_panel, nmds×2, nmds_panel, heatmap ────────
  db<-file.path(output_root,ds,"plots","beta_diversity")
  if(!is.null(beta_plots_all[[ds]])){
    bp<-beta_plots_all[[ds]]
    # PCoA individuais
    if(!is.null(bp$pcoa_plots)) for(dn in names(bp$pcoa_plots)){
      fn<-gsub("[^a-zA-Z0-9]","_",dn)
      tryCatch({save_gg(bp$pcoa_plots[[dn]],file.path(db,paste0("pcoa_",fn)),9,7);ec$plots<-ec$plots+1},
               error=function(e){cat(sprintf("    ERRO pcoa_%s: %s\n",fn,e$message));ec$errors<<-ec$errors+1})}
    # PCoA panel
    if(!is.null(bp$pcoa_panel)){tryCatch({save_gg(bp$pcoa_panel,file.path(db,"pcoa_panel"),16,7);ec$plots<-ec$plots+1},
                                         error=function(e){cat(sprintf("    ERRO pcoa_panel: %s\n",e$message));ec$errors<<-ec$errors+1})}
    # NMDS individuais
    if(!is.null(bp$nmds_plots)&&length(bp$nmds_plots)>0) for(dn in names(bp$nmds_plots)){
      fn<-gsub("[^a-zA-Z0-9]","_",dn)
      tryCatch({save_gg(bp$nmds_plots[[dn]],file.path(db,paste0("nmds_",fn)),9,7);ec$plots<-ec$plots+1},
               error=function(e){cat(sprintf("    ERRO nmds_%s: %s\n",fn,e$message));ec$errors<<-ec$errors+1})}
    # NMDS panel
    if(!is.null(bp$nmds_panel)){tryCatch({save_gg(bp$nmds_panel,file.path(db,"nmds_panel"),16,7);ec$plots<-ec$plots+1},
                                         error=function(e){cat(sprintf("    ERRO nmds_panel: %s\n",e$message));ec$errors<<-ec$errors+1})}
    # Heatmap
    if(!is.null(bp$heatmap)){tryCatch({
      png(file.path(db,"heatmap_top_features.png"),width=14,height=12,units="in",res=300,bg="white")
      draw(bp$heatmap,merge_legend=TRUE); dev.off()
      cat("    Salvo: heatmap_top_features.png\n"); ec$plots<-ec$plots+1
    },error=function(e){tryCatch(dev.off(),error=function(e2)NULL)
      cat(sprintf("    ERRO heatmap: %s\n",e$message)); ec$errors<<-ec$errors+1})}
  }
  
  # ── DA (Fase 11): volcano, barplot, panel ──────────────────────────────────
  dd<-file.path(output_root,ds,"plots","differential_abundance")
  if(!is.null(da_plots_all[[ds]])){
    dap<-da_plots_all[[ds]]
    for(item in list(
      list(p=dap$volcano,n="volcano_maaslin2",    w=11,h=8),
      list(p=dap$bar,    n="barplot_top_da",      w=10,h=9),
      list(p=dap$panel,  n="da_combined_panel",   w=20,h=8)
    )){if(!is.null(item$p)){tryCatch({save_gg(item$p,file.path(dd,item$n),item$w,item$h);ec$plots<-ec$plots+1},
                                     error=function(e){cat(sprintf("    ERRO %s: %s\n",item$n,e$message));ec$errors<<-ec$errors+1})}}
  }
  
  # ── TABELAS EXCEL ──────────────────────────────────────────────────────────
  td<-file.path(output_root,ds,"tables")
  
  # Alpha
  if(!is.null(alpha_stats_all[[ds]])){tryCatch({
    wa<-createWorkbook()
    add_fmt_sheet(wa,"Descriptive",alpha_stats_all[[ds]]$descriptive,"#2E75B6")
    add_fmt_sheet(wa,"Kruskal_Wallis",alpha_stats_all[[ds]]$kruskal_wallis,"#2E75B6")
    add_fmt_sheet(wa,"Wilcoxon",alpha_stats_all[[ds]]$wilcoxon,"#2E75B6")
    add_fmt_sheet(wa,"Per_Sample",alpha_results_all[[ds]],"#2E75B6")
    saveWorkbook(wa,file.path(td,paste0(ds,"_alpha_diversity.xlsx")),overwrite=TRUE)
    cat(paste("    ",ds,"_alpha_diversity.xlsx\n")); ec$tables<-ec$tables+1
  },error=function(e){cat(sprintf("    ERRO alpha xlsx: %s\n",e$message));ec$errors<<-ec$errors+1})}
  
  # Beta
  if(!is.null(beta_stats_all[[ds]])){tryCatch({
    wb<-createWorkbook()
    add_fmt_sheet(wb,"PERMANOVA",beta_stats_all[[ds]]$permanova,"#375623")
    add_fmt_sheet(wb,"ANOSIM",beta_stats_all[[ds]]$anosim,"#375623")
    add_fmt_sheet(wb,"Betadisper",beta_stats_all[[ds]]$betadisper,"#375623")
    add_fmt_sheet(wb,"NMDS_Stress",beta_stats_all[[ds]]$nmds_stress,"#375623")
    saveWorkbook(wb,file.path(td,paste0(ds,"_beta_diversity.xlsx")),overwrite=TRUE)
    cat(paste("    ",ds,"_beta_diversity.xlsx\n")); ec$tables<-ec$tables+1
  },error=function(e){cat(sprintf("    ERRO beta xlsx: %s\n",e$message));ec$errors<<-ec$errors+1})}
  
  # DA
  if(!is.null(da_results_all[[ds]])){tryCatch({
    wd<-createWorkbook()
    add_fmt_sheet(wd,"All_Results",da_results_all[[ds]],"#833C00")
    sig_ds<-da_results_all[[ds]][da_results_all[[ds]]$q_value<q_threshold,]
    if(nrow(sig_ds)>0) add_fmt_sheet(wd,"Significant",sig_ds,"#833C00")
    add_fmt_sheet(wd,"Summary",da_stats_all[[ds]],"#833C00")
    saveWorkbook(wd,file.path(td,paste0(ds,"_differential_abundance.xlsx")),overwrite=TRUE)
    cat(paste("    ",ds,"_differential_abundance.xlsx\n")); ec$tables<-ec$tables+1
  },error=function(e){cat(sprintf("    ERRO da xlsx: %s\n",e$message));ec$errors<<-ec$errors+1})}
  
  cat("\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# 12.2 MASTER CONSOLIDADO
# ─────────────────────────────────────────────────────────────────────────────
cat("=== MASTER CONSOLIDADO ===\n\n")
tryCatch({
  wm<-createWorkbook()
  if(!is.null(all_desc)&&nrow(all_desc)>0) add_fmt_sheet(wm,"Alpha_Descriptive",all_desc,"#2E75B6")
  if(!is.null(all_kw)&&nrow(all_kw)>0) add_fmt_sheet(wm,"Alpha_KW",all_kw,"#2E75B6")
  if(!is.null(all_wx)&&nrow(all_wx)>0) add_fmt_sheet(wm,"Alpha_Wilcoxon",all_wx,"#2E75B6")
  if(!is.null(all_permanova)&&nrow(all_permanova)>0) add_fmt_sheet(wm,"Beta_PERMANOVA",all_permanova,"#375623")
  if(!is.null(all_anosim)&&nrow(all_anosim)>0) add_fmt_sheet(wm,"Beta_ANOSIM",all_anosim,"#375623")
  if(!is.null(all_betadisp)&&nrow(all_betadisp)>0) add_fmt_sheet(wm,"Beta_Betadisper",all_betadisp,"#375623")
  if(!is.null(all_nmds_str)&&nrow(all_nmds_str)>0) add_fmt_sheet(wm,"Beta_NMDS",all_nmds_str,"#375623")
  if(exists("variance_summary")&&!is.null(variance_summary)) add_fmt_sheet(wm,"PCA_Variance",variance_summary,"#7030A0")
  ri2<-do.call(rbind,rarefaction_info_all); if(!is.null(ri2)&&nrow(ri2)>0) add_fmt_sheet(wm,"Rarefaction",ri2,"#7030A0")
  if(!is.null(all_da_stats)&&nrow(all_da_stats)>0) add_fmt_sheet(wm,"DA_Summary",all_da_stats,"#833C00")
  if(!is.null(all_da_sig)&&nrow(all_da_sig)>0) add_fmt_sheet(wm,"DA_Significant",all_da_sig,"#833C00")
  if(!is.null(all_da_results)&&nrow(all_da_results)>0) add_fmt_sheet(wm,"DA_All_Results",all_da_results,"#833C00")
  saveWorkbook(wm,file.path(output_root,"consolidated","MASTER.xlsx"),overwrite=TRUE)
  cat("  MASTER.xlsx salvo!\n"); ec$tables<-ec$tables+1
},error=function(e){cat(sprintf("  ERRO MASTER: %s\n",e$message));ec$errors<-ec$errors+1})

# ─────────────────────────────────────────────────────────────────────────────
# 12.3 SESSION INFO
# ─────────────────────────────────────────────────────────────────────────────
tryCatch({
  sink(file.path(output_root,"session_info.txt"))
  cat(sprintf("Generated: %s\nDir: %s\nGroups: %s\nAges: %s\nDatasets: %s\n\n",
              format(Sys.time()),data_dir,paste(selected_groups,collapse=", "),
              ifelse(is.null(selected_ages),"ALL",paste(selected_ages,collapse=", ")),
              paste(names(datasets_grouped_list),collapse=", ")))
  print(sessionInfo()); sink()
},error=function(e){tryCatch(sink(),error=function(e2)NULL)})

# ─────────────────────────────────────────────────────────────────────────────
# 12.4 AUDITORIA FINAL
# ─────────────────────────────────────────────────────────────────────────────
cat("\n",rep("=",70),"\n  AUDITORIA DE EXPORTACAO\n",rep("=",70),"\n\n")

for(ds in names(datasets_grouped_list)){
  ds_path<-file.path(output_root,ds)
  if(!dir.exists(ds_path)) next
  all_f<-list.files(ds_path,recursive=TRUE)
  npng<-sum(grepl("\\.png$",all_f)); nxlsx<-sum(grepl("\\.xlsx$",all_f))
  cat(sprintf("  [%s] %d PNG + %d XLSX\n",ds,npng,nxlsx))
  for(sub in c("plots/pca","plots/alpha_diversity","plots/beta_diversity",
               "plots/differential_abundance","tables")){
    sp<-file.path(ds_path,sub); if(!dir.exists(sp)) next
    sf<-list.files(sp)
    cat(sprintf("    %s/ (%d)\n",sub,length(sf)))
    for(f in sf) cat(sprintf("      %s\n",f))
  }
  cat("\n")
}

cp<-file.path(output_root,"consolidated")
if(dir.exists(cp)){cf<-list.files(cp);cat(sprintf("  [consolidated] %d\n",length(cf)));for(f in cf) cat(sprintf("    %s\n",f))}
if(file.exists(file.path(output_root,"session_info.txt"))) cat("  session_info.txt\n")

cat(sprintf("\n  TOTAIS: %d plots | %d tabelas | %d erros\n",ec$plots,ec$tables,ec$errors))

# Checklist
cat("\n  === CHECKLIST POR DATASET ===\n")
expected<-c("pca_PC1_vs_PC2","pca_PC1_vs_PC3","sequencing_depth","alpha_all_metrics",
            "alpha_2x2_panel","alpha_observed","alpha_chao1","alpha_shannon","alpha_simpson",
            "pcoa_Bray_Curtis","pcoa_Jaccard","pcoa_panel","nmds_Bray_Curtis","nmds_Jaccard",
            "nmds_panel","heatmap_top_features","volcano_maaslin2","barplot_top_da","da_combined_panel",
            "_alpha_diversity.xlsx","_beta_diversity.xlsx","_differential_abundance.xlsx")

for(ds in names(datasets_grouped_list)){
  cat(sprintf("\n  [%s]\n",ds))
  af<-list.files(file.path(output_root,ds),recursive=TRUE)
  for(exp in expected){
    found<-any(grepl(exp,af,fixed=TRUE))
    cat(sprintf("    %s %s\n",ifelse(found,"OK   ","FALTA"),exp))
  }
}

cat("\n",rep("=",70),"\n  PIPELINE COMPLETO!\n",rep("=",70),"\n")
