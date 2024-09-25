# Helper function 
subset_df <- function(df) {
  return(df[df$p_val_adj < 0.05 & abs(df$avg_log2FC) > 0.3, ])
}

# Run FindMarkers using desired idents/test
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

# Run differential accessibility analysis using Seurat FindMarkers
runDA <- function(seu_obj, test, selected_clusters = NULL, idents) {
  markers_C9HET <- list()
  markers_C9KO <- list()
  DefaultAssay(seu_obj) <- 'peaks'
  Idents(seu_obj) <- idents
  if (!is.null(selected_clusters)){
    seu_obj <- subset(seu_obj, idents = selected_clusters)
  }
  if (is.null(selected_clusters)) {
    selected_clusters = levels(Idents(seu_obj))
  }
  seu_peaks <- DietSeurat(seu_obj, assays = 'peaks') %>%
    SplitObject(split.by = idents)
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

# Annotate DARs to genes using Signac's closest feature
annotateDARs <- function(dar, seu.obj) {
  DefaultAssay(seu.obj) <- 'peaks'
  dar <- lapply(dar, function(genotype) {
    genotype <- lapply(genotype, function(df) {
      if (nrow(df) > 0) {
        df[['genomic_region']] <- rownames(df)
        df <- cbind(df, ClosestFeature(seu.obj, df[['genomic_region']]))
        df$query_region <- NULL # Remove query region column
        return(df)
      }
    })
  })
}
