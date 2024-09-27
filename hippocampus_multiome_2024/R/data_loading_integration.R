@@ -1,377 +0,0 @@
# Sometimes empty fragment objects get created in Seurat object
fragpathfix <- function(oldfrags) {
  t <- 1
  newfrag <- list()
  for (k in seq_along(oldfrags)) {
    if (length(oldfrags[[k]]@cells) != 0){
      newfrag[[t]] <- oldfrags[[k]]
      t <- t + 1
    }
  }
  return(newfrag)
}
# Get genome annotation
genomeAnnotation <- function() {
  annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  seqlevelsStyle(annotation) <- 'UCSC'
  genome(annotation) <- "mm10"
  return(annotation)
}
# Define a function to load RNA data paths, including metadata
load_data <- function() {
  # Temporary variables
  out_list <- list()
  cnames <- c()
  fnames <- c()
  pnames <- c()
  mnames <- c()
  
  # Get all files in the directory
  files.list <- list.files(file_dir, recursive = T)
  
  # Separate files into two categories: counts, metadata
  ## Extract sample names from filenames
  counts <- 
    files.list[grepl("filtered_feature_bc_matrix\\.h5", files.list)]
  fragments <-
    files.list[grepl("atac_fragments\\.tsv\\.gz$", files.list)]
  metadata <-
    files.list[grepl("per_barcode_metrics\\.csv", files.list)]
  peaks <- 
    files.list[grepl("atac_peaks\\.bed", files.list)]
  # Extract sample names from filenames
  for (k in seq_along(fragments)) {
    cnames <- c(cnames, toupper(strsplit(counts[k], '/')[[1]][1]))
    mnames <- c(mnames, toupper(strsplit(metadata[k], '/')[[1]][1]))
    fnames <- c(fnames, toupper(strsplit(fragments[k], '/')[[1]][1]))
    pnames <- c(pnames, toupper(strsplit(peaks[k], '/')[[1]][1]))
  }
  
  # Set names for the lists of files
  names(counts) <- cnames
  names(metadata) <- mnames
  names(fragments) <- fnames
  names(peaks) <- pnames
  
  # Organize files by sample name
  for (sample in mnames) {
    out_list[[sample]] <- list(
      counts = counts[sample],
      metadata = metadata[sample],
      fragments = fragments[sample],
      peaks = peaks[sample]
    )
  }
  return(out_list)
}

prep_data <- function(file.list, file_dir) {
  seurat.list <- list()
  n <- 0
  combined.peaks <- c()
  annotation <- genomeAnnotation()
  for (sample in file.list) {
    n <- n + 1
    # Load counts
    counts <- Read10X_h5(file.path(file_dir, sample$counts))$Peaks
    #counts@Dimnames[[1]] <- gsub(":", "-", counts@Dimnames[[1]])
    # Load peak sets
    peaks <- read.table(
      file = file.path(file_dir, sample$peaks),
      col.names = c("chr", "start", "end")) 
    peaks <- makeGRangesFromDataFrame(peaks)
    # Load metadata
    metadata <- read.table(
      file = file.path(file_dir, sample$metadata),
      stringsAsFactors = FALSE,
      sep = ",",
      header = TRUE,
      row.names = 1
    )[-1, ]
    # Load fragments
    fragments <- CreateFragmentObject(
      path = file.path(file_dir, sample$fragments),
      cells = rownames(metadata),
      validate.fragments = F
    )
    out <- list(
      'counts' = counts,
      'metadata' = metadata,
      'peaks' = peaks,
      'fragments' = fragments
    )
    seurat.list[[names(file.list[n])]] <- out
  }
  # Peak QC, has to be done manually
  combined.peaks <- reduce(x = c(seurat.list$HET1$peaks, seurat.list$HET2$peaks, 
                                 seurat.list$HET3$peaks, seurat.list$HET4$peaks, 
                                 seurat.list$KO1$peaks, seurat.list$KO2$peaks, 
                                 seurat.list$KO3$peaks, seurat.list$KO4$peaks, 
                                 seurat.list$WT1$peaks, seurat.list$WT2$peaks, 
                                 seurat.list$WT3$peaks, seurat.list$WT4$peaks))
  peakwidths <- width(combined.peaks)
  combined.peaks <- combined.peaks[peakwidths < 10000 & peakwidths > 20]
  n <- 0
  for (sample in seurat.list) {
    n <- n + 1
    fmatrix <- FeatureMatrix(
      fragments = sample$fragments,
      features = combined.peaks,
      cells = rownames(sample$metadata)
      )
    cassay <- CreateChromatinAssay(
      counts = fmatrix,
      fragments = sample$fragments,
      annotation = annotation
      )
    seu.obj <- CreateSeuratObject(
      counts = cassay,
      assay = 'ATAC',
      project = 'hippo'
      )
    seurat.list[[cnames[n]]] <- seu.obj[, colnames(sample$counts)]
  }
  return(seurat.list = seurat.list)
}

# Create Seurat objects
add_metadata <- function(seurat.list) {
  new.list <- list()
  n <- 0
  for (seurat.obj in seurat.list) {
    # Set object metadata, except sex which is set manually
    n <- n + 1
    seurat.obj$sample <- names(seurat.list)[[n]]
    seurat.obj$genotype <- paste("C9", sub(".$", '', names(seurat.list)[[n]]), sep = "")
    seurat.obj$region <- "hippocampus"
    seurat.obj <- NucleosomeSignal(seurat.obj)%>%
      TSSEnrichment()
    new.list[[names(seurat.list)[n]]] <- seurat.obj
    # Save Seurat object as an RDS file
    dir_chk('sample_outs')
    saveRDS(seurat.obj,
            file.path("sample_outs", paste(
              names(seurat.list)[n], '_hippo.rds', sep = ""
            )))
  }
  
  new.list[[1]]$sex <- "Male"
  new.list[[2]]$sex <- "Female"
  new.list[[3]]$sex <- "Female"
  new.list[[4]]$sex <- "Male"
  new.list[[5]]$sex <- "Male"
  new.list[[6]]$sex <- "Male"
  new.list[[7]]$sex <- "Female"
  new.list[[8]]$sex <- "Female"
  new.list[[9]]$sex <- "Female"
  new.list[[10]]$sex <- "Male"
  new.list[[11]]$sex <- "Male"
  new.list[[12]]$sex <- "Female"
  return(new.list)
}

# Integration RNA
integrateRNA <- function(merged, num_pcs = 50, resolution = 1.6, k_anchor = 25) {
  # Split by sample and preprocess each subset
  DefaultAssay(merged) <- 'RNA'
  hip.list <- SplitObject(merged, split.by = "sample") %>%
    lapply(FUN = function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
    })
  
  # Select integration features
  features <- SelectIntegrationFeatures(object.list = hip.list)
  
  # Preprocess each subset for integration
  hip.list <- lapply(X = hip.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE, npcs = num_pcs)
  })
  
  # Find integration anchors
  hip.anchors <- FindIntegrationAnchors(object.list = hip.list, 
                                        anchor.features = features,
                                        reduction = "rpca",
                                        k.anchor = k_anchor,
                                        reference = c(1,2,3,4))
  
  # Integrate data
  hip.combined <- IntegrateData(anchorset = hip.anchors, new.assay.name = 'integrated_RNA')
  DefaultAssay(hip.combined) <- "integrated_RNA"
  
  # Additional preprocessing and clustering steps
  hip.combined <- ScaleData(hip.combined, verbose = FALSE)
  hip.combined <- RunPCA(hip.combined, npcs = num_pcs, verbose = FALSE)
  hip.combined <- FindNeighbors(hip.combined, reduction = "pca", dims = 1:num_pcs)
  hip.combined <- RunUMAP(hip.combined, reduction = "pca", dims = 1:num_pcs, umap.method = 'umap-learn', metric = 'correlation')
  hip.combined <- FindClusters(hip.combined, resolution = resolution, algorithm = 4, method = 'igraph', graph.name = 'integrated_RNA_snn')
  return(hip.combined)
}

# Integration peaks
integratePeaks <- function(seu.obj) {
  DefaultAssay(seu.obj) <- 'ATAC'
  # Fix fragments in seurat object, as empty fragments objects duplicate
  # frags <- Fragments(seu.obj) %>%
  #   fragpathfix()
  # Fragments(seu.obj) <- NULL
  # Fragments(seu.obj) <- frags
  atac.obj <- DietSeurat(seu.obj, assay = 'ATAC')
  DefaultAssay(seu.obj) <- 'RNA'
  rna.obj <- DietSeurat(seu.obj, assay = "RNA")
  dir_chk('/mnt/data/single_cell/C9-mice/C9mice_outs/tmp/')
  peaks <- CallPeaks(
    atac.obj,
    group.by = 'class',
    macs2.path = '/home/raghav1881/miniconda3/bin/macs2',
    fragment.tmpdir = '/mnt/data/single_cell/C9-mice/C9mice_outs/tmp',
    outdir = '/mnt/data/single_cell/C9-mice/C9mice_outs/tmp/'
  ) %>%
    keepStandardChromosomes(pruning.mode = "coarse") %>%
    subsetByOverlaps(ranges = blacklist_mm10, invert = TRUE)
  plan(sequential) # Making featurematrix slow on multicore
  macs2_counts <- FeatureMatrix(
    fragments = Fragments(atac.obj),
    features = peaks,
    cells = colnames(atac.obj)
  )
  plan(multisession(workers = 24))
  annotation <- genomeAnnotation()
  atac.obj[["peaks"]] <- CreateChromatinAssay(
    counts = macs2_counts,
    fragments = Fragments(atac.obj),
    annotation = annotation
  )
  
  # Peaks
  DefaultAssay(atac.obj) <- 'peaks'
  peaks.list <-
    SplitObject(DietSeurat(atac.obj, assay = 'peaks'), split.by = "sample")
  peaks.list <- lapply(
    X = peaks.list,
    FUN = function(peaks.obj) {
      DefaultAssay(peaks.obj) <- "peaks"
      peaks.obj <- peaks.obj %>%
        FindTopFeatures(min.cutoff = 5) %>%
        RunTFIDF() %>%
        RunSVD()
    }
  )
  
  # Merge samples
  peaks.combined <- peaks.list[[1]]
  for (i in 2:length(peaks.list)) {
    peaks.combined <- merge(peaks.combined, peaks.list[[i]])
  }
  
  # Process the combined dataset
  peaks.combined <- FindTopFeatures(peaks.combined, min.cutoff = 5)
  peaks.combined <- RunTFIDF(peaks.combined)
  peaks.combined <- RunSVD(peaks.combined)
  peaks.combined <- RunUMAP(peaks.combined, reduction = "lsi", dims = 2:50)
  #p2 <- DimPlot(peaks.combined, group.by = "sample",reduction = 'umap') + ggtitle("ATAC merged")
  
  # compute LSI
  # find integration anchors
  integration.anchors <- FindIntegrationAnchors(
    object.list = peaks.list,
    anchor.features = rownames(peaks.list[[1]]),
    reduction = "rlsi",
    dims = 2:50,
    reference = c(1, 2, 3, 4)
  )
  
  # integrate LSI embeddings
  peaks.integrated <- IntegrateEmbeddings(
    anchorset = integration.anchors,
    reductions = peaks.combined[["lsi"]],
    new.reduction.name = "integrated_lsi",
    dims.to.integrate = 1:50,
  )
  peaks.integrated <-
    RunUMAP(
      peaks.integrated,
      reduction = 'integrated_lsi',
      dims = 2:50,
      reduction.name = "umap.peaks",
      reduction.key = "peaksUMAP_"
    )
  peaks.integrated <-
    FindNeighbors(peaks.integrated, reduction = 'integrated_lsi', dims = 2:50)
  peaks.integrated <-
    FindClusters(
      peaks.integrated,
      reduction = 'integrated_lsi',
      algorithm = 3,
      resolution = 1.1
    )
  p2 <-
    DimPlot_scCustom(
      peaks.integrated,
      reduction = 'umap.peaks',
      group.by = 'seurat_clusters',
      label = T
    ) +
    ggtitle("Peaks integrated") +
    theme(plot.title = element_text(hjust = 0.5)) + NoLegend()
  ggsave('QC-Peaks-Seurat-Clusters.pdf', plot = p2, path = 'qc_figures/')
  p3 <-
    DimPlot_scCustom(
      peaks.integrated,
      reduction = 'umap.peaks',
      group.by = 'class',
      label = T
    ) + 
    ggtitle("Peaks UMAP Annotated") + 
    theme(plot.title = element_text(hjust = 0.5)) + NoLegend()
  ggsave('Peaks_annotated.pdf', plot = p3, path = 'annotated_figures/')
  
  
  integrated_multiome <- hip.combined
  # ATAC related counts should be from 'merge' workflow above
  integrated_multiome[["peaks"]] <- peaks.combined[['peaks']]
  
  # Two reductions: integrated_lsi and umap.peaks are added
  integrated_multiome@reductions <-
    append(
      integrated_multiome@reductions,
      peaks.integrated@reductions,
    )
  integrated_multiome$peaks_clusters <- Idents(peaks.integrated)
  Idents(integrated_multiome) <- 'peaks_clusters'
  saveRDS(integrated_multiome, 'r_objs/integrated_multiome.rds')
  
  return(integrated_multiome)
}

# WNN integration
integrateWNN <- function(integrated_multiome) {
  # WNN
  integrated_multiome <-
    FindMultiModalNeighbors(
      integrated_multiome,
      reduction.list = list("pca", "integrated_lsi"),
      dims.list = list(1:50, 2:50)
    )
  integrated_multiome <- RunUMAP(
    integrated_multiome,
    nn.name = "weighted.nn",
    reduction.name = "umap.wnn",
    reduction.key = "wnnUMAP_",
    return.model = TRUE
  )
  # Supervised projection
  integrated_multiome <-
    RunSPCA(integrated_multiome, assay = 'integrated_RNA', graph = 'wsnn')
  integrated_multiome <- FindClusters(
    integrated_multiome,
    graph.name = "wsnn",
    algorithm = 4,
    method = 'igraph',
    verbose = FALSE,
    resolution = 0.8
  )
  
  return(integrated_multiome)
}
