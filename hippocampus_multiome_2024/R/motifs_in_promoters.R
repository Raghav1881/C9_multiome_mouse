library(Seurat)
library(Signac)
library(cicero)
library(UCell)
set.seed(1234)

get_top_10 <- function(df) {
  df %>% 
    head(10)
}

motifs_in_promoters <- function(obj, motif, deg_list){
  DefaultAssay(obj) <- 'peaks'
  accessible <- Motifs(obj)@data[,motif]
  accessible <- names(accessible)[accessible > 0]
  if (!is.null(accessible)) {  
    accessible_granges <- StringToGRanges(accessible)
    promoters <- promoters(accessible_granges, upstream=2000, downstream=200)
    genes_obj <- Annotation(obj)
    obj_promoters <- promoters(genes_obj, upstream=2000,downstream=200)
    target_genes <- obj_promoters[obj_promoters@ranges %in% promoters@ranges]
    DE_obj <- intersect(rownames(deg_list), target_genes$gene_name)
    if (length(DE_obj) == 0) {return(NULL)}
    return(DE_obj)
  }
}

module_score <- function(obj, features) {
  if (length(features) > 1) {
    mg_scores <- ScoreSignatures_UCell(GetAssayData(obj, slot = 'data'), features = features, name = NULL)
    melted <- reshape2::melt(mg_scores)
    colnames(melted) <- c("Cell","Signature","UCell_score")
    melted$genotype <- obj$genotype
    return(melted)
  } else {
    return(NULL)
  }
}

module_plot <- function(data, clust) {
  genotype_colors <- c("C9WT" = "grey35", "C9HET" = "grey90", "C9KO" = "red")
  
  #Create combined plot
  plt <- ggplot(data, aes(x = Signature, y = UCell_score, group = interaction(Signature, genotype))) +
    geom_violin(aes(fill = genotype), scale = "width", width = 1) +
    geom_boxplot(width = 0.3, position = position_dodge(width = 1), outlier.size = 0) +
    scale_fill_manual(values = genotype_colors) +
    theme_classic() +
    theme(axis.text.x = element_text("Signature"),
          plot.title = element_text(hjust = 0.5)) +
    ggtitle(sprintf(clust))
  
  # plt <- ggviolin(data, x = "Signature", y = "UCell_score", fill = "genotype",
  #                  palette = genotype_colors, width = 1, add = "boxplot", 
  #                 add.params = (list(fill = 'genotype', width = 0.1)))
  return(plt)
}

DefaultAssay(integrated_multiome) <- 'peaks'
clust <- names(deg_filt$markers_C9HET)
degs_C9HET <- deg_filt$markers_C9HET[clust]
degs_C9KO <- deg_filt$markers_C9KO[clust]
motifs_C9HET <- DA_motifs$C9HET[clust]
motifs_C9KO <- DA_motifs$C9KO[clust]
mp_C9HET_filtered <- list()
mp_C9KO <- list()

for (c in clust) {
  sub <- subset(integrated_multiome, idents = c)
  if (nrow(motifs_C9HET[[c]]) != 0) {
    mC9HET <- motifs_C9HET[[c]]$motif
    for (motif in mC9HET) {
      mp_C9HET_filtered[[c]][[motif]] <- motifs_in_promoters(sub, motif, degs_C9HET[[c]])
    }
  }
  if(nrow(motifs_C9KO[[c]]) != 0){
    mC9KO <- motifs_C9KO[[c]]$motif
    for (motif in mC9KO) {
        mp_C9KO[[c]][[motif]] <- motifs_in_promoters(sub, motif, degs_C9KO[[c]])
      }
  }
}

# Remove any clusters that have no DEG overlaps
mp_C9HET_filtered <- Filter(function(x) length(x) > 0, mp_C9HET_filtered)
mp_C9KO <- Filter(function(x) length(x) > 0, mp_C9KO)

mscore_C9HET_RNA <- list()
mscore_C9HET_GAM <- list()
mscore_C9KO_RNA <- list()
mscore_C9KO_GAM <- list()

for (i in seq_along(mp_C9HET_filtered)) {
  DefaultAssay(integrated_multiome) <- 'RNA'
  clust <- names(mp_C9HET_filtered)
  mg_signatures <- mp_C9HET_filtered[[i]]
  mscore_C9HET_RNA[[clust[[i]]]] <- module_score(subset(integrated_multiome,
                                                        idents = clust[[i]]),
                                                 features = mg_signatures)
  DefaultAssay(integrated_multiome) <- 'GAM'
  mscore_C9HET_GAM[[clust[[i]]]] <- module_score(subset(integrated_multiome,
                                                        idents = clust[[i]]),
                                                 features = mg_signatures)
  
}

for (i in seq_along(mp_C9KO)) {
  clust <- names(mp_C9KO)
  mg_signatures <- mp_C9KO[[i]]
  DefaultAssay(integrated_multiome) <- 'RNA'
  mscore_C9KO_RNA[[clust[[i]]]] <- module_score(subset(integrated_multiome,
                                                       idents = clust[[i]]),
                                                features = mg_signatures)
  DefaultAssay(integrated_multiome) <- 'GAM'
  mscore_C9KO_GAM[[clust[[i]]]] <- module_score(subset(integrated_multiome,
                                                       idents = clust[[i]]),
                                                features = mg_signatures)
}
