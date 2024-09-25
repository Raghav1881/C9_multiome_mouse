suppressMessages({
  library(org.Mm.eg.db)
  library(gprofiler2)
})

enrichGOprof <- function(genes, gene_list) {
  # Replace some genes with Ensembl IDs to avoid ambiguous results
  gene_map = c(March1 = "Marchf1",
               Fam155a = "Nalf1",
               Qk = "Qki",
               Atp5o = "Atp5po",
               Hist1h1e = "H1-4",
               Grasp = "Tamalin",
               Lhfp = "Lhfpl6",
               Pak7 = "Pak5",
               Rin2 = "Rassf4",
               Zfp536 = "Znf536",
               Atp5b = "Atp5f1b",
               Atp5e = "Atp5f1e",
               Atp5j = "Atp5pf",
               Atp5h = "Atp5pd",
               Atp5k = "Atp5me",
               Atp5g3 = "Atp5mc3",
               Atp5md = "Atp5mk",
               Atpif1 = "Atp5if1",
               Atp5j2 = "Atp5mf",
               Atp5g1 = "Atp5mc1",
               Zfp804a = "Znf804a",
               Zfp385b = "Znf385b",
               H3f3b = "H3-3a")
  genes <- sapply(genes, function(x) {
    if (x %in% names(gene_map)) {
      gene_map[x]
    } else {
      x
    }
  })
  out <- gost(query = genes, 
              organism = "gp__EywC_IUb0_cwI", 
              custom_bg = gene_list,
              correction_method = "fdr")$result
  if (is.data.frame(out)) {
    out$term_id <-  gsub("%.*", "", out$term_id)
  }
  return(out)
}

# enrich_go <- function(gene) {
#   out <- enrichGO(gene = gene,
#                   universe = gene_list,
#                   OrgDb = org.Mm.eg.db,
#                   ont = "BP",
#                   pAdjustMethod = "BH",
#                   pvalueCutoff  = 0.05,
#                   qvalueCutoff  = 0.1,
#                   readable = T,
#                   keyType = 'SYMBOL')
# }

# simplify_check <- function(out) {
#   if (!is.null(out)) {
#     simplify(out)
#   }
# }

# Organize data to separate upregulated vs downregulated genes #
GOdatalist <- function(seu_obj, de_list) {
  upload_GMT_file(gmtfile = '/home/raghav1881/Downloads/MOUSE_GO_bp_no_GO_iea_symbol.gmt')
  # "gp__EywC_IUb0_cwI" 2024-06-01 MOUSE_GO_bp_no_GO_iea_symbol.gmt
  gene_list <- rownames(seu_obj)
  tmp <- list()
  clustnames_HET <- names(de_list$markers_C9HET)
  clustnames_KO <- names(de_list$markers_C9KO)
  for (k in seq_along(de_list$markers_C9KO)) {
    cc_ko <- clustnames_KO[[k]]
    tmp[[cc_ko]]$ko$upregulated <- subset(de_list[["markers_C9KO"]][[cc_ko]][de_list[["markers_C9KO"]][[cc_ko]]$avg_log2FC > 0.3, ],
                                          select = -c(avg_log2FC)) %>%
      rownames() %>%
      enrichGOprof(gene_list = gene_list)
    tmp[[cc_ko]]$ko$downregulated <- subset(de_list[["markers_C9KO"]][[cc_ko]][de_list[["markers_C9KO"]][[cc_ko]]$avg_log2FC < -0.3, ],
                                            select = -c(avg_log2FC)) %>%
      rownames() %>%
      enrichGOprof(gene_list = gene_list)
  }
  for (k in seq_along(de_list$markers_C9HET)) {
    cc_het <- clustnames_HET[[k]]
    
    tmp[[cc_het]]$het$upregulated <- subset(de_list[["markers_C9HET"]][[cc_het]][de_list[["markers_C9HET"]][[cc_het]]$avg_log2FC > 0.3, ],
                                            select = -c(avg_log2FC)) %>%
      rownames() %>%
      enrichGOprof(gene_list = gene_list)
    tmp[[cc_het]]$het$downregulated <- subset(de_list[["markers_C9HET"]][[cc_het]][de_list[["markers_C9HET"]][[cc_het]]$avg_log2FC < -0.3, ],
                                              select = -c(avg_log2FC)) %>% 
      rownames() %>%
      enrichGOprof(gene_list = gene_list)
  }
  return(tmp)
}

