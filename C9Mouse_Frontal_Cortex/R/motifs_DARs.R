# Function within darMotifs to find accessible motifs
processMotifs <- function(subset.obj, dar_clust, idents_group) {
  top.da.peak <- rownames(dar_clust)
  open.peaks <- AccessiblePeaks(subset.obj, idents = idents_group)
  meta.feature <- GetAssayData(subset.obj, assay = "peaks", slot = "meta.features")
  peaks.matched <- MatchRegionStats(
    meta.feature = meta.feature[open.peaks,],
    query.feature = meta.feature[top.da.peak,],
    n = 50000
  )
  enriched.motifs <- FindMotifs(object = subset.obj,
                                features = top.da.peak)
  return(enriched.motifs)
}

# Add motifs to Seurat object using chromVARmotifs pwms
addMotifs <- function(seu.obj, pfm) {
  DefaultAssay(seu.obj) <- 'peaks'
  seu.obj <- AddMotifs(seu.obj,
                       genome = BSgenome.Mmusculus.UCSC.mm10,
                       pfm = mouse_pwms_v2)
}

# Calculate motif accessibility
darMotifs <- function(seu.obj, DA_lr_filt, idents) {
  Idents(seu.obj) <- idents
  clusts <- levels(seu.obj[[idents]][[idents]])
  out_list <- list()
  for (clust in clusts) {
    subset.obj <- subset(seu.obj, idents = clust)
    Idents(subset.obj) <- 'genotype'
    if (nrow(DA_lr_filt$markers_C9HET[[clust]]) > 0) {
      dar_het <- DA_lr_filt$markers_C9HET[[clust]]
      print(paste("Processing cluster", clust, "for C9HET"))
      out_list$motifs_C9HET[[clust]] <- processMotifs(subset.obj, dar_het, c("C9WT", "C9HET"))
    }
    if (nrow(DA_lr_filt$markers_C9KO[[clust]]) > 0) {
      dar_ko <- DA_lr_filt$markers_C9KO[[clust]]
      print(paste("Processing cluster", clust, "for C9KO"))
      out_list$motifs_C9KO[[clust]] <- processMotifs(subset.obj, dar_het, c("C9WT", "C9KO"))
    }
  }
  return(out_list)
}

da_motif <- function(seu_obj, ident.1) {
  da_motifs <- FindMarkers(
    object = seu_obj,
    ident.1 = ident.1,
    ident.2 = 'C9WT',
    group.by = 'genotype',
    only.pos = T,
    assay = 'chromvar',
    mean.fxn = rowMeans,
    fc.name = 'avg_diff',
  )
  da_motifs$motif <- rownames(da_motifs)
  da_motifs$motif <- gsub('-', '_', da_motifs$motif)
  return(da_motifs)
}

# Function to find motif markers for different subclasses
da_motif_list <- function(seu_obj, selected_clusters, ident) {
  Idents(seu_obj) <- ident
  seu_obj <- subset(seu_obj, idents = selected_clusters)
  DefaultAssay(seu_obj) <- 'chromvar'
  seu.cv <- DietSeurat(seu_obj, assays = c('chromvar'))
  seu.cv[[ident]] <- droplevels(seu.cv[[ident]])
  seu.list <- SplitObject(seu.cv, split.by = ident)
  rm(seu.cv)
  tmp_c9ko <-
    lapply(seu.list, FUN = da_motif, ident.1 = 'C9KO')
  tmp_c9ko <- tmp_c9ko[order(names(tmp_c9ko))]
  tmp_c9het <-
    lapply(seu.list, FUN = da_motif, ident.1 = 'C9HET')
  tmp_c9het <-tmp_c9het[order(names(tmp_c9het))]
  return(list("C9HET" = tmp_c9het, "C9KO" = tmp_c9ko))
}

# Function to annotate motifs according to motif.name
motif_annotation <- function(seu_obj, da_motif_list) {
  motif.names <- Motifs(seu_obj)@motif.names
  
  map_motif_names <- function(motif, motif.names) {
    short_name <- motif.names[[motif]]
    return(short_name)
  }
  
  da_motif_list <- lapply(da_motif_list, function(x) {
    lapply(x, function(y) {
      y$motif.name <- sapply(y$motif, map_motif_names, motif.names = motif.names)
      return(y)
    })
  })
  
  return(da_motif_list)
}
