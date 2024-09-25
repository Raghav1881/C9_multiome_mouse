library(openxlsx)

# Write DEGs master table (all DEGs for all clusters)
mergeDEGs <- function(degs_list, genotype) {
    clustnames <- names(degs_list$markers_C9KO)
    # Add cluster column to each dataframe
    degs_list <- lapply(seq_along(degs_list), function(i) {
        degs_list[[i]]$cluster <- names(degs_list)[[i]]
        return(degs_list[[i]])
    })
    # Merge data frames while retaining the cluster column
    merged_df <- do.call(rbind, degs_list)
    merged_df$genotype <- gsub("markers_", "", genotype)
    return(merged_df)
}

mergeAllDEGs <- function(degs) {
    merged_degs <- lapply(names(degs), function(genotype) {
        if (!is.null(degs[[genotype]])) {
            mergeDEGs(degs[[genotype]], genotype)
        }
    })
    merged_all_degs <- bind_rows(merged_degs)
    return(merged_all_degs)
}

degExcel <- function(deg_merged, degs) {
  wb <- createWorkbook()
  addWorksheet(wb, sheetName = 'master')
  writeDataTable(wb, sheet = 'master', x = deg_merged)
  names(degs$markers_C9HET) <- gsub("/", "_", names(degs$markers_C9HET))
  names(degs$markers_C9KO) <- gsub("/", "_", names(degs$markers_C9KO))
  cnames <- names(degs$markers_C9KO)
  for (clust in cnames) {
    addWorksheet(wb, sheetName = paste(clust, "_C9HET"))
    writeDataTable(wb, sheet = paste(clust, "_C9HET"), x = degs$markers_C9HET[[clust]])
    addWorksheet(wb, sheetName = paste(clust, "_C9KO"))
    writeDataTable(wb, sheet = paste(clust, "_C9KO"), x = degs$markers_C9KO[[clust]])
  }
  saveWorkbook(wb, file = 'DEG/degs.xlsx', overwrite = T)
}

# Create table and check if GO terms exists for genotype/regulation
addTable <- function(nlib, wb, clust, genotype) {
  addWorksheet(wb, sheetName = clust)
  if (!(is.null(nlib[[clust]][[genotype]]$upregulated))) {
    writeData(wb, sheet = clust, x = 'Upregulated')
    writeDataTable(wb, sheet = clust, x = nlib[[clust]][[genotype]]$upregulated@result, startRow = 2)
  }
  if (!(is.null(nlib[[clust]][[genotype]]$downregulated))) {
    writeData(wb, sheet = clust, x = 'Downregulated', startCol = 11)
    writeDataTable(wb, sheet = clust, x = nlib[[clust]][[genotype]]$downregulated@result, startRow = 2, startCol = 11)
  }
}

# Create workbooks for each genotype
GO_excel <- function(nlib) {
  wb_het <- createWorkbook()
  wb_ko <- createWorkbook()
  cnames <- names(nlib)
  for (k in seq_along(nlib)) {
    addTable(nlib, wb = wb_het, clust = cnames[[k]], genotype = 'het')
    addTable(nlib, wb = wb_ko, clust = cnames[[k]], genotype = 'ko')
  }
  # Save workbook for C9HET and C9KO
  saveWorkbook(wb_het, 'GO_analysis/GO_C9HET.xlsx')
  saveWorkbook(wb_ko, 'GO_analysis/GO_C9KO.xlsx')
}
