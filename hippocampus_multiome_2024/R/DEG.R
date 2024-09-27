library(MAST)

# Removes any dataframes with no DEGs (after filtering)
rm_empty_df <- function(df) {
  filtered_df <- df[sapply(df, function(x) nrow(x) > 0)]
  return(filtered_df)
}

# Sort dataframes in names out of order (not used)
sort_df_names <- function(sublist) {
  sorted_names <- names(sublist)[order(names(sublist))]
  return(sublist[sorted_names])
}

# Run MAST using FindMarkers
find_markers_RNA <- function(obj, ident) {
  markers <- FindMarkers(obj,
                         ident.1 = ident,
                         ident.2 = 'C9WT',
                         group.by = 'genotype',
                         min.pct = 0.1,
                         logfc.threshold = 0,
                         test.use = 'MAST',
                         latent.vars = c('sex', 'nCount_RNA'))
  return(markers)
}

# Prepare seurat object to run MAST, error checking in case any issues with idents
runMAST <- function(seu_obj) {
  markers_C9HET <- list()
  markers_C9KO <- list()
  DefaultAssay(seu_obj) <- 'RNA'
  seu_rna <- DietSeurat(seu_obj, assays = 'RNA') %>%
    SplitObject(split.by = 'class')
  ident_list <- names(seu_rna)
  for (l in ident_list) {
    tryCatch({
      markers_C9HET[[l]] <- find_markers_RNA(seu_rna[[l]], 'C9HET')
    }, error = function(e) {
      cat(paste("Error in processing C9HET for class", l, ":", conditionMessage(e), "\n"))
    })
    tryCatch({
      markers_C9KO[[l]] <- find_markers_RNA(seu_rna[[l]], 'C9KO')
    }, error = function(e) {
      cat(paste("Error in processing C9KO for class", l, ":", conditionMessage(e), "\n"))
    })
  }
  markers_C9HET <- rm_empty_df(markers_C9HET)
  markers_C9KO <- rm_empty_df(markers_C9KO)
  out <- list(markers_C9HET = markers_C9HET, markers_C9KO = markers_C9KO) %>%
    lapply(function(sublist) {
      sublist[names(sublist)[order(names(sublist))]]
    })
  return(out)
}

# Find DEGs in Seurat object, return list of DEGs split by genotype
runDEG <- function(hip.combined, idents = 'class', selected_clusters = NULL) {
  Idents(hip.combined) <- idents
  sub <- subset(hip.combined, idents = selected_clusters) %>%
    droplevels()
  out <- runMAST(sub)
  return(out)
}

# Filter DEGs based on log2FC > 0.3 and adjusted p-val < 0.05
filtDEG <- function(degs) {
  degs$markers_C9HET <- lapply(degs$markers_C9HET, function(x){
    x %>%
      filter(abs(avg_log2FC) > 0.3,
             p_val_adj < 0.05) %>%
      mutate(gene = rownames(.))
  })
  degs$markers_C9KO <- lapply(degs$markers_C9KO, function(x){
    x %>%
      filter(abs(avg_log2FC) > 0.3,
             p_val_adj < 0.05) %>%
      mutate(gene = rownames(.))
  })
  return(degs)
}
