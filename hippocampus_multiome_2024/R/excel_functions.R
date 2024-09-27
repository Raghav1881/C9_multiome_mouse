library(openxlsx)

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