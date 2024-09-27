library(Signac)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(chromVARmotifs)
data("mouse_pwms_v2")

# Function to remove empty data frames from a list
rm_empty_df <- function(df) {
  filtered_df <- df[sapply(df, function(x)
    nrow(x) > 0)]
  return(filtered_df)
}

addMotifs <- function(hip.combined, mouse_pwms_v2) {
  # create a peak x motif matrix, where each entry is 1 if the peak contains the motif, otherwise 0
  DefaultAssay(hip.combined) <- 'peaks'
  
  # add motif information
  hip.combined <- AddMotifs(
    object = hip.combined,
    genome = "mm10",
    pfm = mouse_pwms_v2,
    assay = "peaks"
  )
  names(hip.combined@assays$peaks@motifs@pwm) <- hip.combined@assays$peaks@motifs@motif.names
  # run chromvar
  hip.combined <- RunChromVAR(
    object = hip.combined,
    genome = BSgenome.Mmusculus.UCSC.mm10
  )
  return(hip.combined)
}

# Function to find motifs and filter them based on p-value-adjusted
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
  ) %>%
    dplyr::filter(p_val_adj < 0.005, abs(avg_diff) > 0.5)
  da_motifs$motif <- rownames(da_motifs)
  return(da_motifs)
}

# Function to find motif markers for different subclasses
da_motif_list <- function(seu_obj, selected_clusters) {
  Idents(seu_obj) <- "class"
  seu_obj <- subset(seu_obj, idents = selected_clusters)
  DefaultAssay(seu_obj) <- 'chromvar'
  seu.cv <- DietSeurat(seu_obj, assays = c('chromvar'))
  seu.cv$class <- droplevels(seu.cv$class)
  seu.list <- SplitObject(seu.cv, split.by = "class")
  rm(seu.cv)
  tmp_c9ko <-
    lapply(seu.list, FUN = da_motif, ident.1 = 'C9KO')
  tmp_c9ko <- tmp_c9ko[order(names(tmp_c9ko))]
  tmp_c9het <-
    lapply(seu.list, FUN = da_motif, ident.1 = 'C9HET')
  tmp_c9het <-tmp_c9het[order(names(tmp_c9het))]
  return(list("C9HET" = tmp_c9het, "C9KO" = tmp_c9ko))
}
