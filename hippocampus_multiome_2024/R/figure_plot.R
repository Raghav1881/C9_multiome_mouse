library(tidyr)
library(ggrepel)

genotype_colors <- list(C9WT = '#000000',
                        C9HET = '#9436fe',
                        C9KO = '#ff2600')
color_col <- DiscretePalette_scCustomize(num_colors = 28, palette = 'varibow') # Colours for class ident

# Plot for proportions of cells across clusters
prop_plot <- function(hip.combined, ident = NULL, group = NULL) {
  dittoBarPlot(
    object = hip.combined,
    var = ident,
    group.by = group,
    color.panel =  DiscretePalette_scCustomize(num_colors = 28, palette = 'varibow'),
    retain.factor.levels = T)  
}

# Make a plot for cluster markers
markerheatmap <- function(hip.combined, nsforest) {
  markers <- as.character(unlist(nsforest$NSForest_markers))
  col_fun = circlize::colorRamp2(c(-1, 0, 2), viridis(n = 3))
  DefaultAssay(hip.combined) <- 'RNA'
  mat <- as.matrix(AverageExpression(hip.combined, slot = 'data', assays = 'RNA', features = markers, group.by = 'class')$RNA)
  mat <- mat[, order(colnames(mat))]
  clusters <- colnames(mat)
  mat <- apply(mat, 1, scale)
  rownames(mat) <- clusters
  
  ComplexHeatmap::Heatmap(mat, name = "Expression",  
                          cluster_columns = FALSE,
                          cluster_rows = FALSE,
                          show_column_dend = FALSE,
                          show_row_dend = FALSE,
                          col = col_fun,
                          show_column_names = T,
                          use_raster = T,
                          border = T
                          )
  }

# Construct barplots for total number of upregulated and downregulated DEGs   
deg_barplot <- function(deg, selected_clusters) {
  # Merge DEG lists to create one large dataframe
  merge_degs <- function(deg, selected_clusters){
    # Add genotype column in each DEG dataframe, then merge the dataframes
    add_genotype_and_merge <- function(df1, df2, genotype1, genotype2) {
      df1$genotype <- genotype1
      df2$genotype <- genotype2
      merged_df <- dplyr::bind_rows(df1, df2)
      return(merged_df)
    }
    
    # Add a column for the cluster and expression of DEGs
    add_expr_clust <- function(df, cluster_name) {
      df <- df %>%
        mutate(cluster = cluster_name,
              expression = ifelse(avg_log2FC > 0, 'upregulated', 'downregulated'))
      return(df)
    }
    
    merged_data <- list()
    for (cluster in selected_clusters) {
      df1 <- deg$markers_C9HET[[cluster]]
      df2 <- deg$markers_C9KO[[cluster]]
      merged_df <- add_genotype_and_merge(df1, df2, 'C9HET', 'C9KO')
      merged_data[[cluster]] <- merged_df
    }
    merged_data <- merged_data[sapply(merged_data, function(x) !is.null(x))]
    collapsed_data <- bind_rows(lapply(names(merged_data), function(cluster) 
      add_expr_clust(merged_data[[cluster]], cluster)))
    return(collapsed_data)
  }
  collapsed_data <- merge_degs(deg, selected_clusters)
  ### Separate above this line to make dataframe for volcano plots 
  gene_counts <- collapsed_data %>%
    group_by(cluster, genotype, expression) %>%
    summarise(count = n()) %>%
    spread(expression, count, fill = 0)
  summary_counts <- gene_counts %>%
    group_by(cluster, genotype) %>%
    summarise(
      total_upregulated = sum(upregulated),
      total_downregulated = sum(downregulated)
    )
  reshaped_counts <- summary_counts %>%
    pivot_longer(
      cols = c(total_upregulated, total_downregulated),
      names_to = "expression",
      values_to = "count"
    ) %>%
    mutate(regulation = ifelse(
      grepl("upregulated", expression, ignore.case = TRUE),
      "Upregulated",
      "Downregulated"
    ))
  
  plot_upregulated <-
    ggplot(
      reshaped_counts[reshaped_counts$regulation == "Upregulated",],
      aes(
        x = cluster,
        y = count,
        fill = "red",
        color = genotype,
        pattern = genotype
      )
    ) +
    geom_bar_pattern(stat = "identity",
                     position = "dodge",
                     width = 0.7) +
    geom_text(
      aes(label = count),
      position = position_dodge(width = 0.7),
      vjust = -0.5,
      size = 4
    ) +
    labs(title = "Upregulated Genes by Cluster and Genotype",
         y = "Gene Count") +
    theme_prism() +
    theme(axis.title = element_blank(),
          axis.text.x.bottom = element_blank()) +
    scale_shape_manual(values = c("C9HET" = "stripe", "C9KO" = "crosshatch")) +
    scale_color_manual(values = c("C9HET" = "black", "C9KO" = "black"))
  
  plot_downregulated <-
    ggplot(reshaped_counts[reshaped_counts$regulation == "Downregulated",],
           aes(
             x = cluster,
             y = -count,
             fill = 'steelblue2',
             color = genotype,
             pattern = genotype
           )) +
    geom_bar_pattern(stat = "identity",
             position = "dodge",
             width = 0.7,
             fill = 'steelblue2') +
    geom_text(
      aes(label = count),
      position = position_dodge(width = 0.7),
      vjust = 1.5,
      size = 4
    ) +
    labs(x = "Cluster",
         y = "Gene Count") +
    theme_prism() +
    scale_shape_manual(values = c("C9HET" = "stripe", "C9KO" = "crosshatch")) +
    scale_color_manual(values = c("C9HET" = "black", "C9KO" = "black"))
  
  combined_plots <- plot_upregulated / plot_downregulated
  combined_plots + plot_layout(guides = 'auto', axes = "collect_x")
  return(combined_plots)
}

# USED RAW DEG VALUES, use above code to get collapsed dataframe and make volcano plot
volcano_plot <- function(collapsed_data) {
    p <- ggplot(collapsed_data, aes(x = avg_log2FC, y = -log10(p_val_adj), color = genotype)) +
        rasterise(geom_point(aes(color = ifelse(abs(avg_log2FC) <= 0.3 | -log10(p_val_adj) < 1.3, "grey", genotype)),
                             size = 0.75), dpi = 300) +
        labs(title = "Volcano Plot",
             x = "avg_log2FC",
             y = "-log10(p_val_adj)") +
        theme_pubr() +
        facet_wrap(~cluster, scales = "free") +
        scale_color_manual(values = c("#9436fe", "#ff2600", "grey")) +
        geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", color = "black", linewidth = 0.5) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.5)
    
    # Adding labels for the top genes within each cluster facet
    for (clust in unique(collapsed_data$cluster)) {
        subset_data <- filter(collapsed_data, cluster == clust)
        
        top_genes_C9HET <- subset_data %>%
            filter(genotype == "C9HET") %>%
            slice_max(order_by = ifelse(abs(avg_log2FC) > 0.3 & -log10(p_val_adj) > 1.3, -log10(p_val_adj), -Inf), n = 5) %>%
            ungroup()
        
        top_genes_C9KO <- subset_data %>%
            filter(genotype == "C9KO") %>%
            slice_max(order_by = ifelse(abs(avg_log2FC) > 0.3 & -log10(p_val_adj) > 1.3, -log10(p_val_adj), -Inf), n = 5) %>%
            ungroup()
        
        p <- p +
            geom_text_repel(data = top_genes_C9HET,
                            aes(label = gene),
                            box.padding = 0.5,
                            segment.color = 'grey50',
                            force = 2) +
            geom_text_repel(data = top_genes_C9KO,
                            aes(label = gene),
                            box.padding = 0.5,
                            segment.color = 'grey50',
                            force = 2)
    }
    
    return(p)
}

vol_plot <- merge_degs(degs_raw) %>%
  volcano_plot()

# Generate a heatmap of DEGs across samples merging major cluster types
heatmap_DEG <- function(seu, ident, genes) {
  DefaultAssay(seu) <- 'RNA'
  Idents(seu) <- 'celltype'
  seu_sub <- subset(seu, idents = ident) %>%
    ScaleData(features = genes)
  avgexp <- AverageExpression(seu_sub, features = genes, assays = 'RNA', slot = "scale.data", group.by = 'sample', return.seurat = T)
  avgexp$sample <- Idents(avgexp)
}

# Barplot for GO Terms, keep upregulated and downregulated seperately (Enrichment vs depletion), keep fill white and line colour for genotype
GOplot <- function(nlib, colours = NULL, genotype = NULL, celltype = NULL, excit = T) {
  process_nlib <- function(nlib) {
    for (clust in seq_along(nlib)) {
      for (cond_index in seq_along(nlib[[clust]])) {
        for (inner_index in seq_along(nlib[[clust]][[cond_index]])) {
          nlib[[clust]][[cond_index]][[inner_index]]@result$logp10 <- -log10(nlib[[clust]][[cond_index]][[inner_index]]@result$p.adjust)
          nlib[[clust]][[cond_index]][[inner_index]] <- nlib[[clust]][[cond_index]][[inner_index]]@result
        }
      }
    }
    return(nlib)
  }
  combine_dataframes <- function(nlib_excit) {
    add_columns_and_combine <- function(cluster_list, cluster_name) {
      modified_list <- lapply(names(cluster_list), function(treatment_name) {
        lapply(names(cluster_list[[treatment_name]]), function(regulation_name) {
          df <- cluster_list[[treatment_name]][[regulation_name]]
          if (!is.null(df)) {
            df <- df %>%
              mutate(Genotype = treatment_name, Regulation = regulation_name) %>% 
              arrange(p.adjust)
            if (nrow(df) > 0 && any(df$Count > 3)) {
              if (nrow(df) > 0 && any(df$Count > 3)) {
                df <- df[df$Count > 3, ]
                df <- df[1:min(nrow(df), 3), ]
              }
              
            }
          }
        }) %>%
          bind_rows()
      }) %>%
        bind_rows() %>%
        mutate(Cluster = cluster_name)
      return(modified_list)
    }
    combined_dataframes <- lapply(names(nlib_excit), function(cluster_name) {
      cluster_list <- nlib_excit[[cluster_name]]
      combined_cluster <- add_columns_and_combine(cluster_list, cluster_name)
      return(combined_cluster)
    })
    final_combined_dataframe <- bind_rows(combined_dataframes)
    final_split_list <- split(final_combined_dataframe, final_combined_dataframe$Genotype)
    return(final_split_list)
  }
  subset_df <- function(nlib, selected_clusters, regulation){
    nlib_p <- process_nlib(nlib) 
    nlib_p <- subset(nlib_p, names(nlib_p) %in% selected_clusters)
    cdfs <- combine_dataframes(nlib_p)
    collapsed_df <- do.call(rbind, cdfs) %>%
      arrange(
        match(Cluster, selected_clusters),
        match(Genotype, c('het', 'ko')),
        match(Regulation, c('upregulated', 'downregulated')),
        desc(logp10)) %>% 
      subset(Regulation == regulation) %>%
      mutate(row_order = seq_along(logp10))
    #collapsed_df <- mutate(collapsed_df, row_order = seq_along(collapsed_df$logp10))
    return(collapsed_df)
  }
  
  p1 <- ggplot(subsetted_df, aes(x = -row_order,
                              y = logp10,
                              fill = Genotype, 
                              color = Genotype)) +
    geom_bar(stat="identity", position=position_dodge()) +
    geom_text(aes(label = Description),
              hjust = 1.01,
              size = 3,
              color = "black") +
    labs(title = paste("Enrichment of GO Terms in", paste('aa', " C9", toupper(genotype), sep = "")),
         x = "Row Order",
         y = bquote(-log[10]~p[adj])) +
    theme_prism() +
    theme(legend.position = 'bottom',
          axis.text.y = element_blank(),  # Remove y-axis labels
          axis.ticks.y = element_blank(),  # Remove y-axis ticks
          axis.title.y = element_blank(),  # Remove y-axis title
          axis.line.y = element_blank(),  # Set y-axis line color
          axis.line.y.right = element_blank(),  # Remove right y-axis line
          axis.line.x.top = element_blank(),  # Remove top x-axis line
    ) +
    scale_fill_manual(values = alpha(c("het" = "#9436fe", "ko" = '#ff2600'), .2)) +
    scale_color_manual(values = c("het" = "#9436fe", "ko" = '#ff2600')) +
    coord_flip()
  
  p2 <- p1 +
    geom_bar(aes(y = -Count * 1.5 / max(collapsed_df$Count),
                 color = 'grey1'),
             fill = "grey20",
             show.legend = FALSE,
             stat = "identity", position = position_dodge(width = 0.7), alpha = 0.5) +
    geom_text(aes(y = -Count * 1.5 / max(collapsed_df$Count),
                  label = Count),
              position = position_dodge(width = 0.7),
              hjust = 1.5,
              size = 3.5,
              color = "black") +
    theme(axis.text.y = element_blank()) +
    coord_flip()
  
  return(p2)
}
GOplot(collapsed_df, colours = excit_col, genotype = 'het', celltype = 'Excitatory')

# Make heatmaps for DEGs sorted by p-value adjusted
plots_DEGs <- function() {
  sex <- c('Female','Male','Male','Female','Male','Female','Female','Male','Male','Male','Female','Female')
  genotype <- c("WT", "WT", "WT", "WT",
                "C9-HET", "C9-HET", "C9-HET", "C9-HET",
                "C9-KO", "C9-KO", "C9-KO", "C9-KO")
  
  # Order DEGs by p_value_adj
  deg_order <- function(degs_filt, ident, genotype) {
    out <- degs_filt[[genotype]][[ident]] %>%
      arrange(p_val_adj, decreasing = T) %>%
      head(50) %>% 
      pull(gene)
  }
  
  # Generate plots for heatmaps
  generate_heatmap_plot <- function(data, degs, main_title, cluster_cols = F) {
    p <- dittoHeatmap(data[degs, ], show_colnames = TRUE, slot = 'data',
                      main = main_title, border_color = 'black', cluster_cols = cluster_cols,
                      cellheight = 11, cellwidth = 11, annot.by = c('genotype', 'sex'),
                      annot.colors = c('grey30', '#9436fe', '#ff2600', 'steelblue1', 'orange1'), # First 3 colours for WT, C9-HET, C9-KO, then sex
                      heatmap.colors = viridis::viridis(n = 9))
    return(p)
  }
  
  degs_DGhet <- deg_order(degs_filt, 'DG', 'markers_C9HET')
  degs_DGko <- deg_order(degs_filt, 'DG', 'markers_C9KO')
  degs_CA3het <- deg_order(degs_filt, 'CA3', 'markers_C9HET')
  degs_CA3ko <- deg_order(degs_filt, 'CA3', 'markers_C9KO')
  degs_CA1het <- deg_order(degs_filt, 'CA1', 'markers_C9HET')
  degs_CA1ko <- deg_order(degs_filt, 'CA1', 'markers_C9KO')
  
  dg <- subset(hip.combined, idents = 'DG')
  Idents(dg) <- 'sample'
  avgexp_dg <- AverageExpression(dg, 
                                 assays = 'RNA', 
                                 features = c(degs_DGhet, degs_DGko), 
                                 slot = 'data', 
                                 return.seurat = T)
  avgexp_dg$sample <- Idents(avgexp_dg)
  names(genotype)<- levels(avgexp_dg)
  avgexp_dg <- RenameIdents(avgexp_dg, genotype)
  avgexp_dg$genotype <- Idents(avgexp_dg)
  Idents(avgexp_dg) <- 'sample'
  names(sex) <- levels(avgexp_dg)
  avgexp_dg <- RenameIdents(avgexp_dg, sex)
  avgexp_dg$sex <- Idents(avgexp_dg)
  
  ca3 <- subset(hip.combined, idents = 'CA3')
  Idents(ca3) <- 'sample'
  avgexp_ca3 <- AverageExpression(ca3,
                                  assays = 'RNA',
                                  features = c(degs_CA3het, degs_CA3ko),
                                  slot = 'data',
                                  return.seurat = T)
  avgexp_ca3$sample <- Idents(avgexp_ca3)
  names(genotype)<- levels(avgexp_ca3)
  avgexp_ca3 <- RenameIdents(avgexp_ca3, genotype)
  avgexp_ca3$genotype <- Idents(avgexp_ca3)
  Idents(avgexp_ca3) <- 'sample'
  names(sex) <- levels(avgexp_ca3)
  avgexp_ca3 <- RenameIdents(avgexp_ca3, sex)
  avgexp_ca3$sex <- Idents(avgexp_ca3)
  
  ca1 <- subset(hip.combined, idents = 'CA1')
  Idents(ca1) <- 'sample'
  avgexp_ca1 <- AverageExpression(ca1,
                                  assays = 'RNA',
                                  features = c(degs_CA1het, degs_CA1ko),
                                  slot = 'data',
                                  return.seurat = T)
  avgexp_ca1$sample <- Idents(avgexp_ca1)
  names(genotype)<- levels(avgexp_ca1)
  avgexp_ca1 <- RenameIdents(avgexp_ca1, genotype)
  avgexp_ca1$genotype <- Idents(avgexp_ca1)
  Idents(avgexp_ca1) <- 'sample'
  names(sex) <- levels(avgexp_ca1)
  avgexp_ca1 <- RenameIdents(avgexp_ca1, sex)
  avgexp_ca1$sex <- Idents(avgexp_ca1)
  
  # Need to separate male and females for HETs
  Idents(avgexp_dg) <- 'sex'
  avgexp_dgf <- subset(avgexp_dg, idents = 'Female')
  avgexp_dgm <- subset(avgexp_dg, idents = 'Male')
  
  Idents(avgexp_ca3) <- 'sex'
  avgexp_ca3f <- subset(avgexp_ca3, idents = 'Female')
  avgexp_ca3m <- subset(avgexp_ca3, idents = 'Male')
  
  Idents(avgexp_ca1) <- 'sex'
  avgexp_ca1f <- subset(avgexp_ca1, idents = 'Female')
  avgexp_ca1m <- subset(avgexp_ca1, idents = 'Male')
  
  p1 <- generate_heatmap_plot(avgexp_dg, degs_DGhet, 'DG C9-HET')
  p2 <- generate_heatmap_plot(avgexp_dg, degs_DGko, 'DG C9-KO')
  p3 <- generate_heatmap_plot(avgexp_ca3, degs_CA3het, 'CA3 C9-HET')
  p4 <- generate_heatmap_plot(avgexp_ca3, degs_CA3ko, 'CA3 C9-KO')
  p5 <- generate_heatmap_plot(avgexp_ca1, degs_CA1het, 'CA1 C9-HET')
  p6 <- generate_heatmap_plot(avgexp_ca1, degs_CA1ko, 'CA1 C9-KO')
  p7 <- generate_heatmap_plot(avgexp_dgf, degs_DGhet, 'DG C9-HET')
  p8 <- generate_heatmap_plot(avgexp_dgm, degs_DGhet, 'DG C9-HET')
  p9 <- generate_heatmap_plot(avgexp_ca3f, degs_CA3het, 'CA3 C9-HET Female')
  p10 <- generate_heatmap_plot(avgexp_ca3m, degs_CA3het, 'CA3 C9-HET Male')
  p11 <- generate_heatmap_plot(avgexp_ca1f, degs_CA1het, 'CA1 C9-HET Female')
  p12 <- generate_heatmap_plot(avgexp_ca1m, degs_CA1het, 'CA1 C9-HET Male')
}

# CoveragePlot
cov_p <- function(seu_obj, region, ranges, e_upstream, e_downstream, features) {
  sub_range <- granges(seu_obj[ranges, ])
  sub_range$color <- 'orange'
  CoveragePlot(seu_obj, region = region, group.by = 'genotype', features = features,
               extend.upstream = e_upstream, extend.downstream = e_downstream, 
               region.highlight = sub_range,ranges.group.by = 'color') & 
    scale_fill_manual(values = c( "grey35", "#9436fe", "#ff2600"))
  
}

# Plot degs based on logfc values between regions
logFCRegion <- function(degs_1, degs_2, region) {
  degs <- bind_rows(degs_1, degs_2) %>%
    filter(p_val_adj < 0.05) %>%
    filter(region %in% region)
  extract_log2FC <- function(df, region_name) {
    df %>%
      filter(region == region_name) %>%
      select(gene, genotype, avg_log2FC) %>%
      rename_with(~paste0("avg_log2FC_", region_name), avg_log2FC)
  }
  log2FC_hippo <- degs %>% extract_log2FC("hippocampus")
  log2FC_fctx <- degs %>% extract_log2FC("frontal cortex")
  result_df <- full_join(log2FC_hippo, log2FC_fctx, by = c("gene", "genotype"))
  result_df <- result_df %>%
    na.omit()
  return(result_df)
}


logFCGenotype <- function(degs_1, degs_2, genotype1, genotype2) {
  degs <- bind_rows(degs_1, degs_2) %>%
    filter(p_val_adj < 0.05) %>%
    filter(genotype %in% c(genotype1, genotype2))
  
  extract_log2FC <- function(df, genotype_name) {
    df %>%
      filter(genotype == genotype_name) %>%
      select(gene, region, avg_log2FC) %>%
      rename_with(~paste0("avg_log2FC_", genotype_name), avg_log2FC)
  }
  log2FC_hippo <- degs %>% extract_log2FC(genotype1)
  log2FC_fctx <- degs %>% extract_log2FC(genotype2)
  result_df <- full_join(log2FC_hippo, log2FC_fctx, by = c("gene", "region"))
  result_df <- result_df %>%
    na.omit()
  result_df_hippo <- result_df %>% filter(region == "hippocampus")
  result_df_fctx <- result_df %>% filter(region == "frontal cortex")
  return(list(result_df_hippo = result_df_hippo, result_df_fctx = result_df_fctx))
}


plotLogFC <- function(result_df) {
  correlation <- cor(result_df$avg_log2FC_hippocampus, result_df$`avg_log2FC_frontal cortex`)
  
  ggplot(result_df, aes(x = avg_log2FC_hippocampus, y = `avg_log2FC_frontal cortex`, color = genotype)) +
    geom_point(size = 0.3) +
    labs(x = "Avg_log2FC_hippocampus", y = "Avg_log2FC_frontal_cortex", color = "Genotype",
         title = paste("Correlation:", round(correlation, 2))) +
    theme_minimal() + 
    geom_hline(yintercept = 0, color = "darkgrey", size = 0.7, linetype = 'dashed') +
    geom_vline(xintercept = 0, color = "darkgrey", size = 0.7, linetype = 'dashed') +
    facet_wrap(~genotype)
}

# Make heatmaps with this
plots_DEGs <- function(fctx, hip, assay, degs, region, cluster, n_genes) {
  # Helper function to order DEGs by adjusted p-value and select top n genes for a given genotype
  deg_order <- function(degs, n_genes, genotype) {
    out <- degs %>%
      filter(genotype == genotype) %>%
      arrange(p_val_adj) %>%
      head(n = n_genes) %>%
      pull(gene)
    return(out)
  }
  
  # Helper function to generate heatmap plots
  generate_heatmap_plot <- function(data, genes, main_title, cluster_cols = FALSE) {
    p <- dittoHeatmap(data[genes, ], show_colnames = TRUE, slot = 'data',
                      main = main_title, border_color = 'black', cluster_cols = cluster_cols,
                      cellheight = 11, cellwidth = 11, annot.by = c('genotype', 'sex'),
                      annot.colors = c('grey30', '#9436fe', '#ff2600', 'steelblue1', 'orange1'), # WT, C9-HET, C9-KO, then sex
                      heatmap.colors = viridis::viridis(n = 9))
    return(p)
  }
  
  # Define sex and genotype order based on region
  if (region == 'hippocampus') {
    sex <- c('Female', 'Male', 'Male', 'Female', 'Male', 'Female', 'Female', 'Male', 'Male', 'Male', 'Female', 'Female')
    genotype <- c("WT", "WT", "WT", "WT",
                  "C9-HET", "C9-HET", "C9-HET", "C9-HET",
                  "C9-KO", "C9-KO", "C9-KO", "C9-KO")
    region_sub <- subset(hip, idents = cluster)
  } else {
    sex <- c('Male', 'Male', 'Female', 'Female', 'Female', 'Female', 'Female', 'Female', 'Male', 'Male', 'Male', 'Male', 'Female', 'Female')
    genotype <- c("WT", "WT", "WT", "WT", "WT",
                  "C9-HET", "C9-HET", "C9-HET", "C9-HET",
                  "C9-KO", "C9-KO", "C9-KO", "C9-KO")
    region_sub <- subset(fctx, idents = cluster)
  }
  
  # Filter DEGs for the specific region and cluster
  degs_sub <- degs %>% 
    filter(region == region, cluster == cluster)

  # Get top DEGs for C9HET and C9KO
  degs_region_het <- deg_order(degs_sub, n_genes, 'C9-HET')
  degs_region_ko <- deg_order(degs_sub, n_genes, 'C9-KO')
  
  # Calculate average expression for the subset region
  Idents(region_sub) <- 'sample'
  avgexp_region_sub <- AverageExpression(region_sub, 
                                         assays = assay, 
                                         features = c(degs_region_het, degs_region_ko), 
                                         slot = 'data', 
                                         return.seurat = TRUE)
  avgexp_region_sub$sample <- Idents(avgexp_region_sub)
  
  # Rename identities based on genotype and sex
  names(genotype) <- levels(avgexp_region_sub)
  avgexp_region_sub <- RenameIdents(avgexp_region_sub, genotype)
  avgexp_region_sub$genotype <- Idents(avgexp_region_sub)
  
  Idents(avgexp_region_sub) <- 'sample'
  names(sex) <- levels(avgexp_region_sub)
  avgexp_region_sub <- RenameIdents(avgexp_region_sub, sex)
  avgexp_region_sub$sex <- Idents(avgexp_region_sub)
  
  # Generate heatmap plots for C9HET and C9KO
  p1 <- generate_heatmap_plot(avgexp_region_sub, 
                              degs_region_het, 
                              paste(region, '_', cluster, '_', 'C9-HET'))
  p2 <- generate_heatmap_plot(avgexp_region_sub,
                              degs_region_ko,
                              paste(region, '_', cluster, '_', 'C9-KO'))
  
  return(list(C9HET = p1, C9KO = p2))
}



