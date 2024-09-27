# Use cicero to find links across genomic regions 
find_links <- function(subset_obj){
  hip.cds <- as.cell_data_set(x = subset_obj)
  umap_coords <- subset_obj@reductions$umap.peaks
  hip.cicero <- make_cicero_cds(hip.cds, reduced_coordinates = reducedDims(hip.cds)$UMAP)
  genome <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)
  genome.df <- data.frame("chr" = names(genome), "length" = genome)
  conns <- run_cicero(hip.cicero, genomic_coords = genome.df, sample_num = 100)
  
  CCAN_assigns <- generate_ccans(conns)
  links <- ConnectionsToLinks(conns = conns, ccans = CCAN_assigns, threshold = 0.25)
  Links(subset_obj) <- links
  return(subset_obj)                          
}  

oligo <- subset(integrated_multiome)
oligo <- AddMotifs(oligo,
                   genome = BSgenome.Mmusculus.UCSC.mm10,
                   pfm = mouse_pwms_v2)
oligo_dar <- DA_lr_filt$markers_C9KO$Oligo
top.da.peak <- rownames(oligo_dar)
Idents(oligo) <- 'genotype'
open.peaks <- AccessiblePeaks(oligo, idents = c("C9WT", "C9KO"))
meta.feature <- GetAssayData(oligo, assay = "peaks", slot = "meta.features")
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[top.da.peak, ],
  n = 50000
)
enriched.motifs <- FindMotifs(
  object = oligo,
  features = top.da.peak
) %>%
  filter(p.adjust < 0.05)

