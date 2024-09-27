find_links <- function(subset_obj){
  hip.cds <- as.cell_data_set(x = subset_obj)
  umap_coords <- subset_obj@reductions$umap.peaks
  hip.cicero <- make_cicero_cds(hip.cds, reduced_coordinates = reducedDims(hip.cds)
  genome <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)
  genome.df <- data.frame("chr" = names(genome), "length" = genome)
  conns <- run_cicero(hip.cicero, genomic_coords = genome.df, sample_num = 100)
  
  CCAN_assigns <- generate_ccans(conns)
  links <- ConnectionsToLinks(conns = conns, ccans = CCAN_assigns)
  Links(subset_obj) <- links
  return(subset_obj)                          
}  
