suppressMessages({
  library(Seurat)
  library(Signac)
  library(GenomeInfoDb)
  library(EnsDb.Mmusculus.v79)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(dplyr)
  library(scDblFinder)
  library(future)
  library(scCustomize)
  library(reticulate)
  library(BiocParallel)
})

# Set seed and working directory
set.seed(1234)
plan("multicore", workers = 1)
options(future.globals.maxSize = 200000 * 1024^2)
setwd('/media/raghav1881/data_8tb/GitHub/hippocampus_multiome_2024/')
file_dir <- '/mnt/data/single_cell/C9-mice/C9mice_outs/'

# Function to check if directory is not created
dir_chk <- function(path) {
  if (!(dir.exists(path))) {
    dir.create(path)
  }
}

dir_chk('multiome')
dir_chk('r_objs')
setwd('multiome')
source("data_loading_integration.R")
source("DEG.R")
source("GO_analysis.R")
source("DAR.R")

seurat.list <- load_data() %>%
  prep_data() %>%
  add_metadata()

merged <- Merge_Seurat_List(
  list_seurat = seurat.list,
  add.cell.ids = names(seurat.list),
  project = 'C9hippo'
)
merged$sample <- factor(merged$sample, levels = c('WT1','WT2','WT3','WT4','HET1','HET2','HET3','HET4','KO1','KO2','KO3','KO4'))
# Want to retain cells that pass RNA QC and remove cells in one assay but not in other
merged <- merged[, colnames(hip.combined)]
hip.combined <- hip.combined[, colnames(merged)]
hip.combined[['ATAC']] <- merged[['ATAC']]
hip.combined$nCount_ATAC <- merged$nCount_ATAC
hip.combined$nFeature_ATAC <- merged$nFeature_ATAC
hip.combined$nucleosome_signal <- merged$nucleosome_signal
hip.combined$nucleosome_percentile <- merged$nucleosome_percentile
hip.combined$TSS.enrichment <- merged$TSS.enrichment
hip.combined$TSS.percentile <- merged$TSS.percentile

p1 <- VlnPlot_scCustom(
  hip.combined,
  features = c(
    'nCount_ATAC',
    'nFeature_ATAC',
    'TSS.enrichment',
    'nucleosome_signal'
  ),
  group.by = 'sample',
  pt.size = 0,
  raster = F,
  num_columns = 2
)
# Save the plot as a PDF
dir_chk('qc_figures')
ggsave('QC-Sample-Statistics.pdf', plot = p1, path = 'qc_figures/')
saveRDS(merged, 'r_objs/merged.rds')

# QC
hip.combined <- subset(
  x = hip.combined,
  subset = nCount_RNA < 120000 & 
    nCount_RNA > 100 &
    nFeature_RNA > 500 &
    nFeature_RNA < 12000 &
    percent.mt < 5 & 
    nCount_ATAC < 60000 & 
    nCount_ATAC > 500 &
    nFeature_ATAC > 500 &
    nFeature_ATAC < 35000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 2.5 &
    TSS.enrichment < 12 
)

hip.combined <- integratePeaks(hip.combined)
# Remove GNB cluster, as it is extremely noisy in ATAC data
hip.combined <- subset(hip.combined, idents = 'GNB', invert=T)
hip.combined$class <- droplevels(hip.combined$class)

# Reannotate after reintegrating RNA assay
hip.combined <- integrateRNA(hip.combined)
p1 <- DimPlot_scCustom(hip.combined, label = T, group.by = 'seurat_clusters') + theme_prism()
ggsave(filename = 'leiden_clusters.pdf', plot = p1, path = 'qc_figures/')

# Two inhibitory clusters not partitioning, find subcluster to separate
hip.combined <- FindSubCluster(hip.combined, cluster = '21', graph.name = 'integrated_snn', resolution = 0.1, algorithm = 4)
Idents(hip.combined) <- 'sub.cluster'
hip.combined <- FindSubCluster(hip.combined, cluster = '20', graph.name = 'integrated_snn', resolution = 0.1, algorithm = 4)
Idents(hip.combined) <- 'sub.cluster'
hip.combined$sub.cluster <- factor(hip.combined$sub.cluster, levels = str_sort(unique(hip.combined$sub.cluster), numeric = T))
p1 <- DimPlot_scCustom(hip.combined, group.by = 'sub.cluster', label = T) + theme_prism()
ggsave(filename = 'leiden_clusters_subclustered.pdf', plot = p1, path = 'qc_figures/')

# Find Markers
################### PYTHON ########################
runNSForest <- function(hip.combined, cluster_header = 'sub.cluster', n_trees = 100L, n_genes_eval = 6L) {
  use_miniconda("py3.9")
  DefaultAssay(hip.combined) <- "RNA"
  sc <- import("scanpy")
  sp <- import("scipy.sparse")
  sm <- import("scipy.io")
  nf <- import("nsforest")
  
  adata <- sc$AnnData(
    X   = Matrix::t(GetAssayData(hip.combined, slot = 'counts')),
    obs = hip.combined[[]],
    var = GetAssay(hip.combined)[[]]
  )
  adata$obsm$update(umap = Embeddings(hip.combined, "umap"))
  
  NSForest_results <- nf$NSForest(adata, cluster_header = cluster_header, n_trees = n_trees, n_genes_eval = n_genes_eval) 
  dir_chk('NSForest_outputs')
  # Table output between pandas and R data.frame is weird using reticulate
  data.table::fwrite(NSForest_results, 'NSForest_outputs/NSForest_res.csv', quote = F)
  return(NSForest_results)
}
#######################################################
runNSForest(hip.combined, cluster_header = 'sub.cluster')

# Annotate subclusters
sub.clust.ids <- data.table::fread(file = 'NSForest_outputs/NSForest_res.csv')$cluster
names(sub.clust.ids) <- levels(hip.combined$sub.cluster)
Idents(hip.combined) <- "sub.cluster"
hip.combined <- RenameIdents(hip.combined, sub.clust.ids)
hip.combined$subclass <- Idents(hip.combined)
hip.combined$subclass <- factor(Idents(hip.combined), levels = sort(levels(Idents(hip.combined))))
p1 <- DimPlot_scCustom(hip.combined, group.by = 'subclass', pt.size = 0.2, label = T) + theme_prism()
dir_chk('annotated_figures')
ggsave('annotated_subclass.pdf', plot = p1, path = 'annotated_figures/')

# Collapse subclass to class clusters
Idents(hip.combined) <- 'subclass'
hip.combined <- RenameIdents(hip.combined, 'CA1-do' = 'CA1', 'CA1-ve' = 'CA1', 'CA1-ProS' = 'CA1',
                             'CA3-do' = 'CA3', 'CA3-ve' = 'CA3')
hip.combined <- RenameIdents(hip.combined, 'IT HATA' = 'L2/3 IT', 'L2 IT Apr' = 'L2/3 IT', 
                             'L2/3 IT POST-PRE' = 'L2/3 IT', 'L2/3 IT ENTl' = 'L2/3 IT',
                             'L2/3 IT PAR' = 'L2/3 IT', 'L2 IT RSPv POST-PRE' = 'L2/3 IT',
                             'L2 IT Apr' = 'L2/3 IT', 'L3 IT ENTm' = 'L2/3 IT')
hip.combined <- RenameIdents(hip.combined, 'Pvalb Vipr2' = 'Pvalb',
                             'Sst Chodl' = "Sst", 'Sst HPF' = 'Sst')
hip.combined$class <- Idents(hip.combined)
hip.combined$class <- factor(hip.combined$class, levels = sort(levels(hip.combined$class)))
p1 <- DimPlot_scCustom(hip.combined, group.by = 'class', pt.size = 0.2, colors_use = color_col, label = T) + theme_prism()
ggsave('annotated_class.pdf', plot = p1, path = 'annotated_figures/')
#saveRDS(hip.combined, 'r_objs/annotated.rds')
hip.combined <- integratePeaks(hip.combined_new)

# DEG analysis
clusters <- c('Astro-1', 'CA1','CA2','CA3','CA4','DG','Cck','Lamp5','Lamp5 Lhx6','Micro',
              'Oligo','OPC','Pvalb','Sst','Vip')
deg_raw <- runDEG(hip.combined, idents = 'class', selected_clusters = clusters) 
deg_filt <- filtDEG(deg_raw)
# Find DAR with LR and MAST
DA_lr <- runDA(hip.combined, test = 'LR', clusters)
DA_lr_filt <- lapply(DA_lr, function(x) {
  lapply(x, subset_df)
  })
DA_mast <- runDA(hip.combined, test = 'MAST', selected_clusters = clusters)

# Add motif matrix and run chromVAR
hip.combined <- motifAnalysis(hip.combined, mouse_pwms_v2)
DA_motifs <- da_motif_list(integrated_multiome)
dir_chk('motifs')
saveRDS(DA_motifs, 'motifs/DA_motifs.rds')
