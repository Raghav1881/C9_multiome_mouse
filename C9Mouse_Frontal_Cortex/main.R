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
  library(Matrix)
})

# Function to check if directory is not created
dir_chk <- function(path) {
  if (!(dir.exists(path))) {
    dir.create(path)
  }
}

# Sometimes empty fragment objects get created in Seurat object
fragpathfix <- function(obj) {
  current_assay <- obj@active.assay
  DefaultAssay(obj) <- 'ATAC'
  t <- 1
  oldfrags <- Fragments(obj)
  newfrag <- list()
  for (k in seq_along(oldfrags)) {
    if (length(oldfrags[[k]]@cells) != 0) {
      newfrag[[t]] <- oldfrags[[k]]
      t <- t + 1
    }
  }
  Fragments(obj) <- NULL
  Fragments(obj) <- newfrag
  DefaultAssay(obj) <- current_assay
  return(obj)
}

# Remove doublets using scDblFinder
doubletRemoval <- function(seu.obj) {
  sce <- as.SingleCellExperiment(seu.obj, assay = 'RNA')
  sce <- scDblFinder(sce,
                     samples = "sample",
                     clusters = 'seurat_clusters_subsetted',
                     dbr = 0.05)
  seu.obj$scDblFinder.class <- sce$scDblFinder.class
  return(seu.obj)
}

# Find Markers
################### PYTHON ########################
runNSForest <-
  function(seu.obj,
           cluster_header = NULL,
           n_trees = 100L,
           n_genes_eval = 6L) {
    DefaultAssay(seu.obj) <- 'RNA'
    sc <- import("scanpy")
    sp <- import("scipy.sparse")
    sm <- import("scipy.io")
    nf <- import("nsforest")
    
    adata <- sc$AnnData(
      X   = t(GetAssayData(seu.obj, slot = 'counts')),
      obs = seu.obj[[]],
      var = GetAssay(seu.obj)[[]]
    )
    adata$obsm$update(umap = Embeddings(seu.obj, "umap"))
    
    NSForest_results <-
      nf$NSForest(
        adata,
        cluster_header = cluster_header,
        n_trees = n_trees,
        n_genes_eval = n_genes_eval
      )
    dir_chk('NSForest_outputs')
    # Table output between pandas and R data.frame is weird using reticulate
    data.table::fwrite(NSForest_results,
                       'NSForest_outputs/NSForest_res.csv',
                       quote = F)
    return(NSForest_results)
  }
#######################################################

# Get genome annotation
source("R/data_loading_integration.R")
source("R/DEG.R")
source("R/GO_analysis.R")
source("R/DAR.R")
source("R/motifs_DARs.R")
source("R/excel.R")

# Set seed and working directory
set.seed(1234)
plan("multicore", workers = 24)
options(future.globals.maxSize = 20000 * 1024 ^ 2)
setwd('/media/raghav1881/data_8tb/GitHub/C9Mouse_Frontal_Cortex/')
file_dir <- '/media/raghav1881/data_8tb/frontal_cortex/'

dir_chk('r_objs')
seurat.list <- load_data(file_dir)

# Removing unusable samples and reordering
seurat.list <-
  seurat.list[!names(seurat.list) %in% c('HET1', 'HET4')]
seurat.list$HET5 <- seurat.list$HET11
seurat.list$HET6 <- seurat.list$HET12
seurat.list$HET7 <- seurat.list$HET13
seurat.list$HET11 <- NULL
seurat.list$HET12 <- NULL
seurat.list$HET13 <- NULL
seurat.list$WT5 <- seurat.list$`WT4-2`
seurat.list$`WT4-2` <- NULL
sorted_names <- names(seurat.list)[order(names(seurat.list))]
seurat.list <- seurat.list[sorted_names]
seurat.list <- prep_data(seurat.list, file_dir)
seurat.list <- add_metadata(seurat.list)

merged <- Merge_Seurat_List(
  list_seurat = seurat.list,
  add.cell.ids = names(seurat.list),
  project = 'C9ftcx'
)
merged$sample <-
  factor(
    merged$sample,
    levels = c(
      'WT2',
      'WT3',
      'WT4',
      'WT5',
      'WT6',
      'HET2',
      'HET3',
      'HET5',
      'HET6',
      'HET7',
      'KO1',
      'KO2',
      'KO3',
      'KO4'
    )
  )
DefaultAssay(merged) <- 'RNA'
merged[["percent.mt"]] <-
  PercentageFeatureSet(object = merged, pattern = "^mt-")
p1 <- VlnPlot_scCustom(
  merged,
  features = c(
    'nCount_RNA',
    'nFeature_RNA',
    'percent.mt',
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
ggsave(
  'QC-Sample-Statistics.pdf',
  plot = p1,
  path = 'qc_figures/',
  width = 8.5,
  height = 11,
  units = "in"
)
saveRDS(merged, 'r_objs/merged.rds')
merged <- subset(
  x = hip.merged,
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

merged <- doubletRemoval(merged)
p1 <-
  DimPlot_scCustom(merged,
                   label = F,
                   pt.size = 0.2,
                   group.by = 'scDblFinder.class') + theme_prism()
ggsave('predicted_doublets.pdf', plot = p1, path = 'qc_figures/')
Idents(merged) <- 'scDblFinder.class'
merged <- subset(merged, idents = 'singlet') %>%
  droplevels()
merged <- fragpathfix(merged)
DefaultAssay(merged) <- "RNA"
merged <- integrateRNA(merged)
p1 <- DimPlot_scCustom(merged, label = T) + theme_prism()
ggsave(filename = 'leiden_clusters.pdf',
       plot = p1,
       path = 'qc_figures/')

DefaultAssay(merged) <- 'RNA'
runNSForest(merged, cluster_header = 'seurat_clusters')

# Remove cluster 44 as it has no clear markers, very low expressing
p1 <- Cluster_Highlight_Plot(merged, '44') + theme_prism()
dir_chk('qc_figures/noise_clusters')
ggsave('cluster_44_noise.pdf', plot = p1, path = 'qc_figures/noise_clusters/')
p1 <-
  FeaturePlot_scCustom(merged, 'nCount_RNA') + FeaturePlot_scCustom(merged, 'nFeature_RNA') &
  theme_prism()
ggsave(
  'count_features.pdf',
  plot = p1,
  path = 'qc_figures/noise_clusters/',
  width = 8,
  height = 10,
  units = 'in'
)
merged <- subset(merged, idents = '44', invert = T) %>%
  droplevels()
merged <- fragpathfix(merged)
DefaultAssay(merged) <- 'RNA'
merged <- integrateRNA(merged)

# Remove clusters 25 and 26 for thalamic, 34 for both thalamic/noise
dir_chk('qc_figures/thalamic_clusters')
p1 <-
  DimPlot_scCustom(merged, group.by = 'seurat_clusters', label = T) + theme_prism()
ggsave(filename = 'leiden_clusters.pdf',
       plot = p1,
       path = 'qc_figures/thalamic_clusters/')
p1 <-
  FeaturePlot_scCustom(
    merged,
    c('Drd1', 'Ido1','Ppp1r1b','Adora2a','Gpr83','Penk','Drd1','Dclk3','Mhrt')
  )
ggsave(
  filename = 'thalamic_clusters.pdf',
  plot = p1,
  path = 'qc_figures/thalamic_clusters/',
  width = 12,
  height = 12
)
saveRDS(merged, 'r_objs/with_thalamic.rds')
merged <-
  subset(merged, idents = c('25', '26', '34'), invert = T) %>%
  droplevels() %>%
  fragpathfix()
DefaultAssay(merged) <- 'RNA'
merged <- integrateRNA(merged)

# Find subclusters for 23, 37, 35
merged <-
  FindSubCluster(
    merged,
    cluster = '23',
    graph.name = 'integrated_RNA_snn',
    resolution = 0.1,
    algorithm = 4
  )
Idents(merged) <- 'sub.cluster'
merged <-
  FindSubCluster(
    merged,
    cluster = '37',
    graph.name = 'integrated_RNA_snn',
    resolution = 0.1,
    algorithm = 4
  )
Idents(merged) <- 'sub.cluster'
merged <-
  FindSubCluster(
    merged,
    cluster = '35',
    graph.name = 'integrated_RNA_snn',
    resolution = 0.1,
    algorithm = 4
  )
merged$sub.cluster <-
  factor(merged$sub.cluster, levels = str_sort(unique(merged$sub.cluster), numeric = T))
p1 <-
  DimPlot_scCustom(merged,
                   group.by = 'sub.cluster',
                   label = T,
                   pt.size = 0.1) + theme_prism()
ggsave(filename = 'leiden_clusters_subclustered_with37_3.pdf',
       plot = p1,
       path = 'qc_figures/')
p1 <- Cluster_Highlight_Plot(merged, '37_3', pt.size = 0.1)
ggsave(filename = '37_3noise.pdf',
       plot = p1,
       path = 'qc_figures/')

# Remove cluster 37_1 as it is from LSX (Ano1, Myo5b) and 37_3 as there are no discernable markers and it does not partition
merged <- subset(merged, idents = c('37_3', '37_1'), invert = T) %>%
  droplevels() %>%
  fragpathfix()
DefaultAssay(merged) <- 'RNA'
merged <- integrateRNA(merged)
p1 <-
  DimPlot_scCustom(merged, group.by = 'seurat_clusters', label = T) + theme_prism()
ggsave(filename = 'leiden_clusters.pdf',
       plot = p1,
       path = 'qc_figures/')

#Find subcluster for 35, 38 and 42
merged <-
  FindSubCluster(
    merged,
    cluster = '35',
    graph.name = 'integrated_RNA_snn',
    resolution = 0.2,
    algorithm = 4
  )
Idents(merged) <- 'sub.cluster'
merged <-
  FindSubCluster(
    merged,
    cluster = '38',
    graph.name = 'integrated_RNA_snn',
    resolution = 0.1,
    algorithm = 4
  )
Idents(merged) <- 'sub.cluster'
merged <-
  FindSubCluster(
    merged,
    cluster = '42',
    graph.name = 'integrated_RNA_snn',
    resolution = 0.1,
    algorithm = 4
  )
Idents(merged) <- 'sub.cluster'
merged$sub.cluster <-
  factor(merged$sub.cluster, levels = str_sort(unique(merged$sub.cluster), numeric = T))
p1 <-
  DimPlot_scCustom(merged, group.by = 'sub.cluster', label = T) + theme_prism()
ggsave(filename = 'leiden_subclusters_final.pdf',
       plot = p1,
       path = 'qc_figures/')
saveRDS(merged, 'r_objs/merged_without_thalamic_LHX.rds')
runNSForest(merged, cluster_header = 'sub.cluster')

# Annotate sub clusters
anno <- fread('NSForest_outputs/NSForest_res.csv')$cluster
Idents(merged) <- 'sub.cluster'
names(anno) <- levels(merged$sub.cluster)
merged <- RenameIdents(merged, anno)
merged$subclass <- Idents(merged)
merged$subclass <-
  factor(merged$subclass, levels = sort(levels(merged$subclass)))

# Remove OECs as they are olfactory bulb cells
Idents(merged) <- 'subclass'
merged <- subset(merged, idents = 'OEC', invert = T) %>%
  droplevels()
color_col <-
  DiscretePalette_scCustomize(
    num_colors = 25,
    palette = 'varibow',
    seed = 1234,
    shuffle_pal = T
  ) # Colours for plot
p1 <-
  DimPlot_scCustom(
    merged,
    group.by = 'subclass',
    pt.size = 0.2,
    colors_use = color_col,
    label = T
  ) + theme_prism()
dir_chk('annotated_figures')
ggsave('annotated_subclass.pdf', plot = p1, path = 'annotated_figures/')

# Find DEGs
deg_raw <- runDEG(merged, idents = 'subclass')
deg_filt <- filtDEG(deg_raw)
dir_chk('DEG')
deg_merged <- mergeAllDEGs(deg_filt)
degExcel(deg_merged, deg_filt)
saveRDS(deg_raw, 'DEG/deg_raw.rds')
saveRDS(deg_filt, 'DEG/deg_filt.rds')

# Run GO analysis
nlib <- GOdatalist(merged, deg_filt)
# Filter go terms with fewer than 3 genes
nlib_filt <- filtGO(nlib)

dir_chk('GO_analysis')
GO_excel(nlib)
saveRDS(nlib, 'GO_analysis/nlib.rds')

# Merge peaks
merged <- integratePeaks(merged)
# DA analysis
DA_lr <- runDA(merged, test = 'LR', idents = 'subclass')
DA_lr_filt <- lapply(DA_lr, function(x) {
  lapply(x, subset_df)
})
dir_chk('DAR')
saveRDS(DA_lr, 'DAR/DA_lr.rds')
saveRDS(DA_lr_filt, 'DAR/DA_lr_filt.rds')

# Add motifs to Seurat object
merged <- addMotifs(merged, pfm = mouse_pwms_v2)
# Find accessible motifs in DARs
motif_outs_raw <-
  darMotifs(merged, DA_lr_filt = DA_lr_filt, idents = 'subclass')
motif_outs_filt <-
  lapply(motif_outs_raw, function(x) {
    lapply(x, function(y) {
      y <- filter(y, pvalue < 0.05)
    })
  })
dir_chk("motif")
saveRDS(motif_outs_filt, 'motif/motif_outs_filt.rds')
# Run chromvar to calculate motif accessibility
merged <- RunChromVAR(object = merged,
                      genome = BSgenome.Mmusculus.UCSC.mm10)
saveRDS(merged, 'r_objs/merged_multiome_motifs.rds')

# Run chromvar
merged <- RunChromVAR(object = merged,
                      genome = BSgenome.Mmusculus.UCSC.mm10)
da_motifs_chromvar <- da_motif_list(merged,
                                    selected_clusters = levels(merged$subclass),
                                    ident = 'subclass') %>%
  motif_annotation(seu_obj = merged)
da_motifs_chromvar_filt <- lapply(da_motifs_chromvar, function(x) {
  lapply(x, function(y) {
    y_filtered <- y[y$p_val_adj < 0.05, ]
    return(y_filtered)
  })
})

saveRDS(da_motifs_chromvar, 'motif/motifs_chromvar.rds')
saveRDS(da_motifs_chromvar_filt, 'motif/motifs_chromvar_filt.rds')
