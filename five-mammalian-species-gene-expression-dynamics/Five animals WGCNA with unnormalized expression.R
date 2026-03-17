library(WGCNA)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(igraph)
library(ggraph)
library(clusterProfiler)
library(org.Hs.eg.db)
allowWGCNAThreads()  # Enable parallel processing for WGCNA

#Read expr 
setwd("F:/Single cell analysis/Parallel development/Human/Objects")
Human_stemcell_regress <- read_rds("Human_stemcell_regress_25_04_29.rds")

setwd("F:/Single cell analysis/Parallel development/Macaque/Objects")
Macaque_stemcell_regress <- read_rds("Macaque_stemcell_regress_25_05_19.rds")

setwd("F:/Single cell analysis/Parallel development/Ferret/Objects")
Ferret_stemcell_regress <- read_rds("Ferret_stemcell_regress_25_05_19.rds")

setwd("F:/Single cell analysis/Parallel development/Mouse/Objects")
Mouse_stemcell_regress <- read_rds("Mouse_stemcell_regress_25_05_19.rds")

setwd("F:/Single cell analysis/Parallel development/Paper figures/Figure 7/objects")
Gene_atlas <- read.csv("Atlas conserved age genes with linear model and HAR and disease.csv", header = T)

setwd("F:/Single cell analysis/Mouse and rat/Network and trend/objects")
Rat_stemcell_regress <- read_rds("Rat_stemcell_regress_25_08_06.rds")

GOI <- c("AXIN2", "PAX6","HES1","SOX2","WNT3","WNT4","WNT5A","WNT5B","WNT7A","WNT7B", "WNT8B", "FRZB", "WLS", "LMO2")
GOI <- c("HMGA2", "CCND1", "IGF2BP1","IGF2BP2","HMGA1", "CRYAB", "PDK1")

ligand_gene <- c("WNT1", "WNT2", "WNT2B", "WNT3A", "WNT4", "WNT5A", "WNT5B", "WNT6","WNT7B", "WNT8A", "WNT8B", "WNT9A", "WNT9B", "WNT10A", "WNT10B", "WNT11", "WNT16", "PORCN", "WLS", "PGAP1", "GPC3")
receptor_gene <- c("FZD1", "FZD2", "FZD3", "FZD4", "FZD5", "FZD6", "FZD8", "FZD9", "FZD10", "LRP5", "LRP6", "RYK", "ROR2", "PTK7", 'PKD1', 'TSPAN12')
GOI <- ligand_gene
GOI <- receptor_gene
GOI <- sort(c(ligand_gene, receptor_gene))

GOI <- unique(GO_retrieve("GO:0045880"))     #SHH positive
GOI <- unique(c(GO_retrieve("GO:0030177")))  #Wnt positive
GOI <- unique(c(GO_retrieve("GO:0032008")))  #mTOR positive regulation
GOI <- unique(c(GO_retrieve("GO:0045747")))  #NOTCH positive
GOI <- unique(c(GO_retrieve("GO:0045743")))  #FGF positive
GOI <- unique(c(GO_retrieve("GO:0008543")))  #FGF 
GOI <- unique(c(GO_retrieve("GO:0030513")))  #BMP positive
FindVariableFeatures()

# Helper function to compute average expression by cell type and stage
get_avg_expr_by_celltype_stage <- function(seurat_obj, cell_types, genes) {
  # Get genes present in the SCT assay
  sct_genes <- rownames(seurat_obj[["SCT"]])
  genes_present <- intersect(genes, sct_genes)
  
  if (length(genes_present) == 0) {
    stop("None of the provided genes are found in the SCT assay.")
  }
  
  if (length(genes_present) < length(genes)) {
    missing_genes <- setdiff(genes, genes_present)
    warning(paste("The following genes were not found in the SCT assay and will be set to zero:", paste(missing_genes, collapse = ", ")))
  }
  
  # Compute average expression for present genes
  avg_expr <- AverageExpression(
    seurat_obj,
    assays = "SCT",
    layer = "data",  # For Seurat v5
    features = genes_present,
    group.by = c("stemcell_type", "stage_order")
  )$SCT
  
  cell_counts <- table(seurat_obj$stemcell_type, seurat_obj$stage_order)
  # Convert dgMatrix to matrix
  avg_mat <- as.matrix(avg_expr)
  
  cell_counts_df <- as.data.frame(cell_counts)
  colnames(cell_counts_df) <- c("stemcell_type", "stage_order", "cell_count")
  cell_counts_df$group_id <- paste(cell_counts_df$stemcell_type, cell_counts_df$stage_order, sep = "_")
  low_count_groups <- cell_counts_df$group_id[cell_counts_df$cell_count < 50]
  
  cols_to_zero <- colnames(avg_mat)[colnames(avg_mat) %in% low_count_groups]
  avg_mat[, cols_to_zero] <- 0
  
  # Get all unique stages and cell types
  all_stages <- unique(seurat_obj$stage_order)
  all_cell_types <- cell_types
  
  # Generate all possible combinations of cell types and stages
  all_combinations <- expand.grid(cell_type = all_cell_types, stage = all_stages)
  all_col_names <- paste(all_combinations$cell_type, all_combinations$stage, sep = "_")
  
  # Initialize a matrix with zeros for all possible combinations and all input genes
  full_avg_mat <- matrix(0, nrow = length(genes), ncol = length(all_col_names))
  rownames(full_avg_mat) <- genes
  colnames(full_avg_mat) <- all_col_names
  
  # Map the computed average expression to the corresponding columns and genes
  computed_col_names <- colnames(avg_mat)
  for (col in computed_col_names) {
    # Extract cell type and stage from the column name
    parts <- strsplit(col, "_")[[1]]
    ct <- parts[1]
    stage <- parts[2]
    target_col <- paste(ct, stage, sep = "_")
    if (target_col %in% all_col_names) {
      # Assign only the present genes
      full_avg_mat[genes_present, target_col] <- avg_mat[, col]
    }
  }
  
  # Verify the expected number of columns
  expected_cols <- length(all_stages) * length(cell_types)
  if (ncol(full_avg_mat) != expected_cols) {
    warning("Some stage-cell type combinations may be missing in the final matrix.")
  }
  
  return(full_avg_mat)
}

# Function to retrieve expression matrix across species
Retrive_expr_mat <- function(genes) {
  cell_types <- c("vRG", "ExN-IPC", "oRG", "tRG", "Glia-IPC", "Astro-imma", "OPC")
  species <- c("Human", "Macaque", "Ferret", "Mouse", "Rat")
  
  # Subset gene atlas
  gene_atlas_subset <- Gene_atlas[Gene_atlas$Conserved.genes %in% genes, ]
  
  # Extract species-specific genes
  genes_human <- gene_atlas_subset$Human.gene.symbol
  genes_macaque <- gene_atlas_subset$Macaque.gene.symbol
  genes_ferret <- gene_atlas_subset$Ferret.gene.symbol
  genes_mouse <- gene_atlas_subset$Mouse.gene.symbol
  genes_rat <- gene_atlas_subset$Mouse.gene.symbol
  
  # Compute average expression for each species
  human_avg <- get_avg_expr_by_celltype_stage(Human_stemcell_regress, cell_types, genes_human)
  macaque_avg <- get_avg_expr_by_celltype_stage(Macaque_stemcell_regress, cell_types, genes_macaque)
  ferret_avg <- get_avg_expr_by_celltype_stage(Ferret_stemcell_regress, cell_types, genes_ferret)
  mouse_avg <- get_avg_expr_by_celltype_stage(Mouse_stemcell_regress, cell_types, genes_mouse)
  rat_avg <- get_avg_expr_by_celltype_stage(Rat_stemcell_regress, cell_types, genes_rat)
  
  #rownames(ferret_avg)[rownames(ferret_avg) == "WNT3"] <- "WNT3A"
  # Adjust gene names for consistency
  adjusted_rownames_human <- rownames(human_avg)
  adjusted_rownames_macaque <- rownames(macaque_avg)
  adjusted_rownames_ferret <- sub("\\.BLAST.*", "", rownames(ferret_avg))
  adjusted_rownames_mouse <- toupper(rownames(mouse_avg))
  adjusted_rownames_rat <- toupper(rownames(rat_avg))
  
  # Find common genes across species
  common_adjusted_rownames <- Reduce(intersect, list(
    adjusted_rownames_human,
    adjusted_rownames_macaque,
    adjusted_rownames_ferret,
    adjusted_rownames_mouse,
    adjusted_rownames_rat
  ))
  
  # Subset to common genes
  human_avg <- human_avg[adjusted_rownames_human %in% common_adjusted_rownames, , drop = FALSE]
  macaque_avg <- macaque_avg[adjusted_rownames_macaque %in% common_adjusted_rownames, , drop = FALSE]
  ferret_avg <- ferret_avg[sub("\\.BLAST.*", "", rownames(ferret_avg)) %in% common_adjusted_rownames, , drop = FALSE]
  mouse_avg <- mouse_avg[toupper(rownames(mouse_avg)) %in% common_adjusted_rownames, , drop = FALSE]
  rat_avg <- rat_avg[toupper(rownames(rat_avg)) %in% common_adjusted_rownames, , drop = FALSE]
  
  # Set consistent row names
  rownames(human_avg) <- common_adjusted_rownames
  rownames(macaque_avg) <- common_adjusted_rownames
  rownames(ferret_avg) <- common_adjusted_rownames
  rownames(mouse_avg) <- common_adjusted_rownames
  rownames(rat_avg) <- common_adjusted_rownames
  
  # Add species prefix to columns
  colnames(human_avg) <- paste("human", colnames(human_avg), sep = "_")
  colnames(macaque_avg) <- paste("macaque", colnames(macaque_avg), sep = "_")
  colnames(ferret_avg) <- paste("ferret", colnames(ferret_avg), sep = "_")
  colnames(mouse_avg) <- paste("mouse", colnames(mouse_avg), sep = "_")
  colnames(rat_avg) <- paste("rat", colnames(rat_avg), sep = "_")
  
  # Combine into a single expression matrix
  expr_mat <- cbind(human_avg, macaque_avg, ferret_avg, mouse_avg, rat_avg)
  return(expr_mat)
}

# WGCNA Analysis Setup
GOI_expr <- Retrive_expr_mat(GOI)  # Replace GOI with your ~1000 genes
datExpr <- as.matrix(t(GOI_expr))  # Samples as rows, genes as columns
rownames(datExpr) <- colnames(GOI_expr)  # Ensure proper row names
colnames(datExpr)

# Define traits
samples <- rownames(datExpr)
split_samples <- strsplit(samples, "_")
species_list <- sapply(split_samples, `[`, 1)
cell_type_list <- sapply(split_samples, `[`, 2)
stage_list <- as.numeric(sapply(split_samples, `[`, 3))

datTraits <- data.frame(
  species = species_list,
  cell_type = cell_type_list,
  stage = stage_list,
  row.names = samples  # Explicitly set row names
)

# Choose soft-thresholding power
powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
power <- sft$powerEstimate
if (is.na(power)) power <- 6  # Default power if no good fit
#power <- 6  # 6 is actually pretty good


# Construct network and identify modules
net <- blockwiseModules(
  datExpr,
  power = power,
  TOMType = "unsigned",
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  verbose = 3
)

# Get module colors and eigengenes
#moduleColors <- labels2colors(net$colors)
MEs <- moduleEigengenes(datExpr, net$colors)$eigengenes
colnames(MEs)

#plot_me_trends("grey", MEs, datTraits, cell_type_colors)

# Create species-cell type specific stage traits
unique_combos <- unique(datTraits[, c("species", "cell_type")])
trait_cols <- list()
for (i in 1:nrow(unique_combos)) {
  sp <- unique_combos$species[i]
  ct <- unique_combos$cell_type[i]
  trait_name <- paste(sp, ct, "stage", sep = "_")
  trait_vals <- ifelse(datTraits$species == sp & datTraits$cell_type == ct, datTraits$stage, NA)
  trait_cols[[trait_name]] <- trait_vals
}
datTraits <- cbind(datTraits, do.call(cbind, trait_cols))

# Correlate module eigengenes with traits
moduleTraitCor <- cor(MEs, datTraits[, -(1:3)], use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Identify modules with opposite trends
cell_types <- c("vRG", "ExN-IPC", "oRG", "tRG", "Glia-IPC", "Astro-imma", "OPC")
modules_of_interest <- list()
for (ct in cell_types) {
  human_trait <- paste("human", ct, "stage", sep = "_")
  other_traits <- paste(c("macaque", "ferret", "mouse", "rat"), ct, "stage", sep = "_")
  
  mods <- which(
    moduleTraitCor[, human_trait] < -0.3 &
      apply(moduleTraitCor[, other_traits], 1, min) > 0.3 &
      moduleTraitPvalue[, human_trait] < 0.05 &
      apply(moduleTraitPvalue[, other_traits], 1, max) < 0.05
  )
  
  if (length(mods) > 0) {
    modules_of_interest[[ct]] <- names(MEs)[mods]
  }
}


modules_of_interest

module <- "yellow"
module <- "brown"
module <- "turquoise"
module <- "green"
module <- "blue"
module <- "grey"
names(net$colors)[net$colors == module]
length(names(net$colors)[net$colors == module])

intersect(names(net$colors)[net$colors == module], tf_list)
write.csv(names(net$colors)[net$colors == "blue"], "module blue.csv")
write.csv(names(net$colors)[net$colors == "turquoise"], "module turquoise.csv")
write.csv(names(net$colors)[net$colors == "brown"], "module brown.csv")
write.csv(names(net$colors)[net$colors == "grey"], "module grey.csv")


###############################################################################
# Plot module eigengene trends
plot_me_trends <- function(module, MEs, datTraits, cell_type_colors, stages_to_plot = NULL) {
  me_col <- paste0("ME", module)
  plot_data <- data.frame(
    ME = MEs[, me_col],
    species = datTraits$species,
    cell_type = datTraits$cell_type,
    stage = datTraits$stage
  )
  
  # Reorder species factor levels
  plot_data$species <- factor(plot_data$species, levels = c("human", "macaque", "mouse", "rat", "ferret"))
  
  # Filter plot_data based on stages_to_plot if provided
  if (!is.null(stages_to_plot)) {
    # Note: stages_to_plot should be a list (e.g., list(human = 1:8, macaque = 1:4)),
    # and the stages must match the type of datTraits$stage (numeric, factor, etc.)
    plot_data_list <- split(plot_data, plot_data$species)
    plot_data_filtered <- lapply(names(plot_data_list), function(sp) {
      if (sp %in% names(stages_to_plot)) {
        stages <- stages_to_plot[[sp]]
        plot_data_list[[sp]] %>% dplyr::filter(stage %in% stages)
      } else {
        plot_data_list[[sp]]  # Keep all stages if species not in stages_to_plot
      }
    })
    plot_data <- do.call(rbind, plot_data_filtered)
  }
  
  # Calculate y-axis limits based on the data range
  y_limits <- range(plot_data$ME, na.rm = TRUE)
  
  p <- ggplot(plot_data, aes(x = stage, y = ME, color = cell_type, group = cell_type)) +
    geom_smooth(aes(fill = cell_type), method = "loess", se = TRUE, alpha = 0.05, size = 0.8, span = 2, level = 0.5) +
    facet_wrap(~ species, scales = "free_x", ncol = 4) +
    scale_color_manual(values = cell_type_colors) +
    scale_fill_manual(values = cell_type_colors) +
    coord_cartesian(ylim = y_limits) +
    labs(
      title = paste("Module Eigengene Trends for Module", module),
      subtitle = "Across Species and Cell Types",
      x = "Developmental Stage",
      y = "Module Eigengene Expression",
      color = "Cell Type",
      fill = "Cell Type"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# Usage example
cell_type_colors <- c("vRG" = "#ddc49e", "ExN-IPC" = "#e98074", "oRG" = "#46a166", "tRG" = "#82d83d", 
                      "Glia-IPC" = "#a6c4d1", "Astro-imma" = "#1e7894", "OPC" = "#897dcb")

stages_to_plot <- list(human = 1:17, macaque = 1:8, ferret = 1:8, mouse = 1:11, rat = 1:3)
#stages_to_plot <- list(human = 1:17, macaque = 4:8, ferret = 1:8, mouse = 1:11)
plot_me_trends("turquoise", MEs, datTraits, cell_type_colors, stages_to_plot)
plot_me_trends("yellow", MEs, datTraits, cell_type_colors)
plot_me_trends("brown", MEs, datTraits, cell_type_colors)
plot_me_trends("green", MEs, datTraits, cell_type_colors)
plot_me_trends("blue", MEs, datTraits, cell_type_colors, stages_to_plot)
plot_me_trends("grey", MEs, datTraits, cell_type_colors)




###############################################################################
# Plot co-expression network with integrated core TF identification
plot_coexpression_network <- function(module, datExpr, power, tf_list, top_connections = 50, label_top_n = 10) {
  module_genes <- names(net$colors)[net$colors == module]
  TOM <- TOMsimilarityFromExpr(datExpr[, module_genes], power = power)
  
  # Select top connections
  TOM_upper <- TOM[upper.tri(TOM)]
  threshold <- quantile(TOM_upper, 1 - (top_connections / choose(length(module_genes), 2)))
  links <- which(TOM > threshold, arr.ind = TRUE)
  links <- links[links[,1] < links[,2], ]
  links_df <- data.frame(from = module_genes[links[,1]], to = module_genes[links[,2]])
  
  if (nrow(links_df) == 0) {
    stop("No connections above the threshold. Try increasing top_connections.")
  }
  
  # Create the network graph
  network <- graph_from_data_frame(links_df, directed = FALSE)
  
  # Get the genes actually present in the network
  network_genes <- V(network)$name
  
  # Calculate module eigengenes
  colors <- rep(module, length(module_genes))
  MEs <- moduleEigengenes(datExpr[, module_genes], colors)$eigengenes
  MM <- cor(datExpr[, module_genes], MEs, use = "p")
  
  # Subset MM to only the genes in the network
  MM_network <- MM[network_genes, , drop = FALSE]
  
  # Identify all TFs in the module that are also in the network
  module_tfs <- intersect(module_genes, tf_list)
  tfs_in_network <- module_tfs[module_tfs %in% network_genes]
  
  # Add node attributes
  V(network)$kME <- MM_network[, 1]
  V(network)$size <- 2 + 10 * (V(network)$kME - min(V(network)$kME)) / (max(V(network)$kME) - min(V(network)$kME))
  V(network)$color <- scales::col_numeric("Blues", domain = V(network)$kME)(V(network)$kME)
  V(network)$label <- NA
  
  # Highlight and label all TFs in the network
  if (length(tfs_in_network) > 0) {
    tf_indices <- which(V(network)$name %in% tfs_in_network)
    V(network)$color[tf_indices] <- "red"
    V(network)$label[tf_indices] <- V(network)$name[tf_indices]
  }
  
  # Label top 10 genes by kME that are in the network (excluding TFs if desired)
  top_genes <- module_genes[order(MM[, 1], decreasing = TRUE)]
  top_genes_in_network <- top_genes[top_genes %in% network_genes & !(top_genes %in% tfs_in_network)][1:label_top_n]
  top_idx <- which(V(network)$name %in% top_genes_in_network)
  V(network)$label[top_idx] <- V(network)$name[top_idx]
  
  # Plot with ggraph
  p <- ggraph(network, layout = "fr") +
    geom_edge_link(alpha = 0.4, color = "grey50", edge_width = 0.5) +
    geom_node_point(aes(size = size, color = color)) +
    geom_node_text(aes(label = label), repel = TRUE, size = 3, fontface = "bold", color = "black") +
    scale_color_identity() +
    scale_size_identity() +
    theme_void() +
    labs(
      title = paste("Co-expression Network for Module", module),
      subtitle = paste("TFs in red | Top", top_connections, "Connections"),
      caption = "Node size and color intensity reflect kME (module membership)"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      plot.caption = element_text(hjust = 0.5, size = 10, color = "grey50"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  return(p)
}


# Usage example
tf_list <- read.csv("Human_TF.csv")$Symbol
plot_coexpression_network("yellow", datExpr, power, tf_list, top_connections = 500)
plot_coexpression_network("brown", datExpr, power, tf_list, top_connections = 200, label_top_n = 5)
plot_coexpression_network("turquoise", datExpr, power, tf_list, top_connections = 200, label_top_n = 5)
plot_coexpression_network("green", datExpr, power, tf_list, top_connections = 200)
plot_coexpression_network("blue", datExpr, power, tf_list, top_connections = 200, label_top_n = 5)
plot_coexpression_network("grey", datExpr, power, tf_list, top_connections = 200, label_top_n = 5)


###############################################################################
plot_whole_coexpression_network <- function(datExpr, power, net, tf_list, gene_sets, top_connections = 50) {
  library(ggraph)
  library(igraph)
  library(scales)
  
  # Get all module genes
  all_genes <- colnames(datExpr)
  TOM <- TOMsimilarityFromExpr(datExpr, power = power)
  
  # Select top connections across all genes
  TOM_upper <- TOM[upper.tri(TOM)]
  threshold <- quantile(TOM_upper, 1 - (top_connections / choose(length(all_genes), 2)))
  links <- which(TOM > threshold, arr.ind = TRUE)
  links <- links[links[,1] < links[,2], ]
  links_df <- data.frame(from = all_genes[links[,1]], to = all_genes[links[,2]])
  
  if (nrow(links_df) == 0) {
    stop("No connections above the threshold. Try increasing top_connections.")
  }
  
  # Create the network graph
  network <- graph_from_data_frame(links_df, directed = FALSE)
  network_genes <- V(network)$name
  # Assign module colors to nodes
  V(network)$module <- net$colors[match(V(network)$name, names(net$colors))]
  
  # Identify core TFs for each module
  core_tfs <- list()
  for (mod in unique(net$colors)) {
    module_genes <- names(net$colors)[net$colors == mod]
    module_tfs <- intersect(module_genes, tf_list)
    tfs_in_network <- module_tfs[module_tfs %in% network_genes]
    if (length(module_tfs) > 0) {
      MEs <- moduleEigengenes(datExpr[, module_genes], rep(mod, length(module_genes)))$eigengenes
      MM <- cor(datExpr[, module_genes], MEs, use = "p")
      tf_MM <- MM[module_tfs, , drop = FALSE]
      core_tfs[[mod]] <- rownames(tf_MM)[which.max(tf_MM[,1])]
    }
  }
  
  
  # Add node attributes
  V(network)$color <- V(network)$module
  V(network)$label <- NA
  V(network)$border <- "transparent"
  
  # Highlight core TFs
  for (mod in names(core_tfs)) {
    core_tf <- core_tfs[[mod]]
    if (core_tf %in% V(network)$name) {
      core_idx <- which(V(network)$name == core_tf)
      V(network)$color[core_idx] <- "red"
      V(network)$label[core_idx] <- core_tf
    }
  }
  module_tfs <- intersect(module_genes, tf_list)
  tfs_in_network <- module_tfs[module_tfs %in% network_genes]
  if (length(tfs_in_network) > 0) {
    tf_indices <- which(V(network)$name %in% tfs_in_network)
    V(network)$color[tf_indices] <- "red"
    V(network)$label[tf_indices] <- V(network)$name[tf_indices]
  }
  # Annotate gene sets with colored borders
  for (gene_set in gene_sets) {
    genes <- gene_set[[1]]
    color <- gene_set[[2]]
    set_idx <- which(V(network)$name %in% genes)
    V(network)$border[set_idx] <- color
    V(network)$label[set_idx] <- V(network)$name[set_idx]
  }
  
  # Layout to group modules
  layout <- create_layout(network, layout = "fr", weights = E(network)$weight)
  
  # Plot with ggraph
  p <- ggraph(layout) +
    geom_edge_link(alpha = 0.4, color = "grey50", edge_width = 0.5) +
    geom_node_point(aes(color = color), size = 5, stroke = 1.5, shape = 21, fill = V(network)$color, color = V(network)$border) +
    geom_node_text(aes(label = label), repel = TRUE, size = 3, fontface = "bold", color = "black") +
    scale_color_identity() +
    theme_void() +
    labs(
      title = "Whole Co-expression Network",
      subtitle = "Modules and Core TFs",
      caption = "Core TFs in red; gene sets with colored borders"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      plot.caption = element_text(hjust = 0.5, size = 10, color = "grey50"),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  return(p)
}

RPL_gene <- unique(c(c("RPL36","RPL36AL","RPL37A","RPL21","RPL17","RPL7","RPL39","RPL32","RPL29","RPL22L1"), GO_retrieve("GO:0022625")))
NDUF_gene <- c("NDUFS5","NDUFA6","NDUFA3","NDUFA2","NDUFB4","NDUFB2","COX6C", "COX6B1", "COX5B","NDUFA13")
CILIUM_gene <- c("CFAP46","CFAP52","CFAP126","DRC1","WDR63","DAW1")
gene_sets <- list(
  list(CILIUM_gene, "pink")
)

plot_whole_coexpression_network(datExpr, power, net, tf_list, gene_sets, top_connections = 2000)


###############################################################################
# Plot dendrogram
plot_dendrogram <- function(net) {
  plotDendroAndColors(
    net$dendrograms[[1]], 
    net$colors[net$blockGenes[[1]]],
    "Module Colors",
    dendroLabels = F, 
    hang = 0.2,
    addGuide = TRUE, 
    guideHang = 0.05,
    main = "Gene Dendrogram and Module Assignments",
    cex.main = 1.2,
    font.main = 2
  )
}

plot_dendrogram(net)



###############################################################################
# Plot module-trait heatmap
plot_module_trait_heatmap <- function(moduleTraitCor, moduleTraitPvalue, MEs, trait_cols) {
  textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) <- dim(moduleTraitCor)
  
  par(mar = c(6, 8.5, 3, 3))
  labeledHeatmap(
    Matrix = moduleTraitCor, 
    xLabels = trait_cols,
    yLabels = names(MEs),
    colors = blueWhiteRed(50),
    textMatrix = textMatrix, 
    setStdMargins = FALSE,
    cex.text = 0.6, 
    zlim = c(-1, 1),
    main = "Module-Trait Relationships",
    cex.main = 1.2,
    font.main = 2
  )
}

# Usage example
trait_cols <- grep("_stage$", names(datTraits), value = TRUE)
plot_module_trait_heatmap(moduleTraitCor, moduleTraitPvalue, MEs, trait_cols)



###############################################################################
#Plot heatmap directly from df
wgcna_heatmap <- function(datExpr, genes, cell_types) {
  # Parse row names to extract species, cell type, and stage
  genes <- intersect(colnames(datExpr), genes)
  row_info <- data.frame(
    row_name = rownames(datExpr),
    stringsAsFactors = FALSE
  )
  row_info$species <- sapply(strsplit(row_info$row_name, "_"), `[`, 1)
  row_info$cell_type <- sapply(strsplit(row_info$row_name, "_"), `[`, 2)
  row_info$stage <- as.integer(sapply(strsplit(row_info$row_name, "_"), `[`, 3))
  
  # Filter rows based on specified cell types
  filtered_row_info <- row_info[row_info$cell_type %in% cell_types, ]
  
  # Define species order and set factor levels
  species_order <- c("rat", "mouse", "human", "macaque",  "ferret")
  filtered_row_info$species <- factor(filtered_row_info$species, levels = species_order)
  filtered_row_info$cell_type <- factor(filtered_row_info$cell_type, levels = cell_types)
  
  # Sort samples by species, cell type, and stage
  sorted_row_info <- filtered_row_info[order(
    filtered_row_info$species,
    filtered_row_info$cell_type,
    filtered_row_info$stage
  ), ]
  
  # Get sorted row names
  sorted_rows <- sorted_row_info$row_name
  
  # Subset and transpose datExpr: rows as genes, columns as samples
  heatmap_data <- t(datExpr[sorted_rows, genes, drop = FALSE])
  
  # Create column split factor for species and cell type
  split_factor <- paste(sorted_row_info$species, sorted_row_info$cell_type, sep = "_")
  split_factor <- factor(split_factor, levels = unique(split_factor))
  
  # Define color palette from white to dark green
  heatmap_colors <- colorRampPalette(c("white", "darkgreen"))(50)
  
  # Create annotation for the top of the heatmap
  ha <- HeatmapAnnotation(
    Group = split_factor,
    col = list(Group = setNames(rainbow(length(unique(split_factor))), unique(split_factor))),
    show_legend = TRUE,
    show_annotation_name = FALSE
  )
  
  # Create the heatmap with annotation
  ht <- Heatmap(
    heatmap_data,
    name = "Expression",
    col = heatmap_colors,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    column_split = split_factor,
    column_gap = unit(2, "mm"),
    border = FALSE,
    rect_gp = gpar(col = NA),  # Remove cell borders
    top_annotation = ha  # Add the annotation bar
  )
  
  return(ht)
}

cell_types <- c("vRG")
cell_types <- c("ExN-IPC")
cell_types <- c("Astro-imma")
cell_types <- c("oRG")
cell_types <- c("OPC")

ht <- wgcna_heatmap(datExpr, setdiff(GOI, c("HMGA2","CCND1", "WNT7A", "LMO2","WNT3")), cell_types)
ht <- wgcna_heatmap(datExpr,c("HMGA2"), cell_types)
ht <- wgcna_heatmap(datExpr,c("CCND1"), cell_types)
ht <- wgcna_heatmap(datExpr,GOI, cell_types)
draw(ht)


##################################################################
#simply plot expression trend of specific gene
plot_df_trends <- function(datExpr, gene_to_plot, datTraits, cell_type_colors, stages_to_plot = NULL) {
  # Check if the gene exists in datExpr
  if (!gene_to_plot %in% colnames(datExpr)) {
    stop("Gene not found in datExpr.")
  }
  
  # Extract expression data for the specified gene
  expr_data <- datExpr[, gene_to_plot, drop = FALSE]
  
  # Combine with traits data
  plot_data <- data.frame(
    expression = expr_data[, 1],
    species = datTraits$species,
    cell_type = datTraits$cell_type,
    stage = datTraits$stage
  )
  
  # Reorder species factor levels
  plot_data$species <- factor(plot_data$species, levels = c("rat", "mouse", "ferret","macaque","human"))
  
  # Filter plot_data based on stages_to_plot if provided
  if (!is.null(stages_to_plot)) {
    plot_data_list <- split(plot_data, plot_data$species)
    plot_data_filtered <- lapply(names(plot_data_list), function(sp) {
      if (sp %in% names(stages_to_plot)) {
        stages <- stages_to_plot[[sp]]
        plot_data_list[[sp]] %>% dplyr::filter(stage %in% stages)
      } else {
        plot_data_list[[sp]]  # Keep all stages if species not in stages_to_plot
      }
    })
    plot_data <- do.call(rbind, plot_data_filtered)
  }
  
  # Calculate y-axis limits based on the data range
  y_limits <- range(plot_data$expression, na.rm = TRUE)
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = stage, y = expression, color = cell_type, group = cell_type)) +
    geom_smooth(aes(fill = cell_type), method = "loess", se = TRUE, alpha = 0.05, size = 0.8, span = 1.5, level = 0.5) +
    facet_wrap(~ species, scales = "free_x", ncol = 5) +
    scale_color_manual(values = cell_type_colors) +
    scale_fill_manual(values = cell_type_colors) +
    coord_cartesian(ylim = y_limits) +
    labs(
      title = paste("Expression Trends for Gene", gene_to_plot),
      subtitle = "Across Species and Cell Types",
      x = "Developmental Stage",
      y = "Expression Level",
      color = "Cell Type",
      fill = "Cell Type"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# Usage example
stages_to_plot <- list(human = 1:15, macaque = 1:15, ferret = 1:15, mouse = 1:15, rat = 1:3)
stages_to_plot <- list(human = 1:5, macaque = 2:4, ferret = 1:3, mouse = 3:7, rat = 1:3)
#cell_type_colors <- c("vRG" = "#ddc49e", "ExN-IPC" = "red", "oRG" = "white", "tRG" = "white", "Glia-IPC" = "white", "Astro-imma" = "white", "OPC" = "white")
cell_type_colors <- c("vRG" = "#ddc49e", "ExN-IPC" = "#e98074", "oRG" = "#46a166", "tRG" = "#82d83d", 
                      "Glia-IPC" = "#a6c4d1", "Astro-imma" = "#1e7894", "OPC" = "#897dcb")
#GOI <- c("AXIN2", "WNT3","WNT4","WNT5A","WNT5B","WNT7A","WNT7B", "WNT8B", "FRZB", "WLS", "LMO2", "HMGA2", "CCND1")
plot_df_trends(datExpr, "ST18", datTraits, cell_type_colors, stages_to_plot)
plot_df_trends(datExpr, "AXIN2", datTraits, cell_type_colors, stages_to_plot)
plot_df_trends(datExpr, "PAX6", datTraits, cell_type_colors, stages_to_plot)
plot_df_trends(datExpr, "SOX2", datTraits, cell_type_colors, stages_to_plot)
plot_df_trends(datExpr, "IGF2BP1", datTraits, cell_type_colors, stages_to_plot)
plot_df_trends(datExpr, "HMGA2", datTraits, cell_type_colors, stages_to_plot)

plot_df_trends(datExpr, "WNT7B", datTraits, cell_type_colors, stages_to_plot)
plot_df_trends(datExpr, "WNT7A", datTraits, cell_type_colors, stages_to_plot)
plot_df_trends(datExpr, "WNT3", datTraits, cell_type_colors, stages_to_plot)
plot_df_trends(datExpr, "WNT4", datTraits, cell_type_colors, stages_to_plot)
plot_df_trends(datExpr, "WLS", datTraits, cell_type_colors, stages_to_plot)
plot_df_trends(datExpr, "FRZB", datTraits, cell_type_colors, stages_to_plot)
