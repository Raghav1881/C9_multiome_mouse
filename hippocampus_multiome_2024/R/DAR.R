subset_df <- function(df) {
  return(df[df$p_val_adj < 0.05 & abs(df$avg_log2FC) > 0.3, ])
}

find_markers <- function(obj, ident, test) {
  markers <- FindMarkers(obj,
                         ident.1 = ident,
                         ident.2 = 'C9WT',
                         group.by = 'genotype',
                         min.pct = 0.01,
                         logfc.threshold = 0,
                         test.use = test,
                         latent.vars = c('sex', 'nCount_peaks'))
  return(markers)
}

runDA <- function(seu_obj, test, selected_clusters) {
  markers_C9HET <- list()
  markers_C9KO <- list()
  DefaultAssay(seu_obj) <- 'peaks'
  Idents(seu_obj) <- 'class'
  seu_obj <- subset(seu_obj, idents = selected_clusters)
  seu_peaks <- DietSeurat(seu_obj, assays = 'peaks') %>%
    SplitObject(split.by = 'class')
  ident_list <- names(seu_peaks)
  for (l in ident_list) {
    tryCatch({
      markers_C9HET[[l]] <- find_markers(seu_peaks[[l]], 'C9HET', test)
    }, error = function(e) {
      cat(paste("Error in processing C9HET for class", l, ":", conditionMessage(e), "\n"))
    })
    tryCatch({
      markers_C9KO[[l]] <- find_markers(seu_peaks[[l]], 'C9KO', test)
    }, error = function(e) {
      cat(paste("Error in processing C9KO for class", l, ":", conditionMessage(e), "\n"))
    })
  }
  out <- list(markers_C9HET = markers_C9HET, markers_C9KO = markers_C9KO) %>%
    lapply(function(sublist) {
      sublist[names(sublist)[order(names(sublist))]]
    })
  return(out)
}
