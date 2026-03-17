library(WGCNA)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(igraph)
library(ggraph)
library(clusterProfiler)
library(org.Hs.eg.db)

setwd("F:/Single cell analysis/Mouse and rat/Network and trend/objects")

human_mean_expr <- read_rds("Human_vRG_expr_mean.rds") 
macaque_mean_expr <- read_rds("Macaque_vRG_expr_mean.rds") 
mouse_mean_expr <- read_rds("Mouse_vRG_expr_mean.rds") 
ferret_mean_expr <- read_rds("Ferret_vRG_expr_mean.rds") 
rat_mean_expr <- read_rds("Rat_vRG_expr_mean.rds") 

human_mean_expr <- read_rds("Human_ExN_IPC_expr_mean.rds") 
macaque_mean_expr <- read_rds("Macaque_ExN_IPC_expr_mean.rds") 
mouse_mean_expr <- read_rds("Mouse_ExN_IPC_expr_mean.rds") 
ferret_mean_expr <- read_rds("Ferret_ExN_IPC_expr_mean.rds") 
rat_mean_expr <- read_rds("Rat_ExN_IPC_expr_mean.rds") 

GOI <- c("AXIN2", "WNT3","WNT4","WNT5A","WNT5B","WNT7A","WNT7B", "WNT8B", "FRZB", "WLS", "LMO2")
GOI <- c("HMGA2","CCND1", "HES1","NOTCH2")
GOI <- GO_retrieve("GO:0030177") #positive regulation of WNT
intersect(sort(read.csv("module_v2blue.csv")$x), GO_retrieve("GO:0016055"))
GOI <- read.csv("module_blue.csv")$x
GOI <- read.csv("module_brown.csv")$x
ligand_gene <- c("WNT1", "WNT2", "WNT2B", "WNT3", "WNT3A", "WNT4", "WNT5A", "WNT5B", "WNT6", "WNT7A", "WNT7B", "WNT8A", "WNT8B", "WNT9A", "WNT9B", "WNT10A", "WNT10B", "WNT11", "WNT16", "PORCN", "WLS", "PGAP1", "GPC3")
ligand_gene <- c( "WNT5A", "WNT5B", "WNT7B", "WLS", "PGAP1", "GPC3")
receptor_gene <- c("FZD1", "FZD3", "FZD4", "FZD5", "FZD6", "FZD8", "FZD9", "FZD10", "LRP5", "LRP6",  "ROR2", "PTK7", 'PKD1')
GOI <- ligand_gene
GOI <- receptor_gene

GOI <- unique(GO_retrieve("GO:0045880"))     #SHH positive
GOI <- unique(c(GO_retrieve("GO:0030177")))  #Wnt positive
GOI <- unique(c(GO_retrieve("GO:0032008")))  #mTOR positive regulation
GOI <- unique(c(GO_retrieve("GO:0045747")))  #NOTCH positive
GOI <- unique(c(GO_retrieve("GO:0045743")))  #FGF positive
GOI <- unique(c(GO_retrieve("GO:0008543")))  #FGF 
GOI <- unique(c(GO_retrieve("GO:0030513")))  #BMP positive

# Helper function to interpolate expression to a common scaled grid
# This function scales stages from 0 to 1 and interpolates gene expression onto a fixed number of bins.
get_scaled_binned_expr <- function(mean_expr, species, num_bins = 30) {
  # Convert rownames to numeric stages and sort them
  stages <- as.numeric(rownames(mean_expr))
  sorted_idx <- order(stages)
  stages <- stages[sorted_idx]
  mean_expr <- mean_expr[sorted_idx, , drop = FALSE]
  
  # Compute min and max stages for scaling
  min_stage <- min(stages)
  max_stage <- max(stages)
  scaled_stages <- (stages - min_stage) / (max_stage - min_stage)
  
  # Create a regular grid from 0 to 1 with num_bins points
  grid <- seq(0, 1, length.out = num_bins)
  
  # Interpolate each gene's expression onto the grid using linear approximation
  expr_grid <- apply(mean_expr, 2, function(y) {
    approx(scaled_stages, y[sorted_idx], xout = grid, method = "linear")$y
  })
  
  # Replace any NA values (from extrapolation) with 0
  expr_grid[is.na(expr_grid)] <- 0
  
  # Transpose: rows = genes, columns = bins; add species prefix to column names
  expr_grid <- t(expr_grid)
  colnames(expr_grid) <- paste(species, seq_len(num_bins), sep = "_")
  
  return(expr_grid)
}

# Function to retrieve expression matrix across species (expanded to five species, no cell types)
# This prepares a combined matrix of interpolated expressions for genes of interest.
# Handles genes that may not be present in all species by setting missing expressions to 0.
Retrive_expr_mat <- function(genes = NULL, num_bins = 30) {
  species <- c("human", "macaque", "mouse", "ferret", "rat")
  mean_expr_list <- list(
    human = human_mean_expr,
    macaque = macaque_mean_expr,
    mouse = mouse_mean_expr,
    ferret = ferret_mean_expr,
    rat = rat_mean_expr
  )
  
  # Determine all genes to include (union if genes=NULL, or provided genes)
  if (is.null(genes)) {
    all_genes <- unique(unlist(lapply(mean_expr_list, colnames)))
  } else {
    all_genes <- genes
    # Warn if some genes are not found in any species
    found_genes <- unique(unlist(lapply(mean_expr_list, colnames)))
    missing <- setdiff(genes, found_genes)
    if (length(missing) > 0) {
      warning(paste("Genes not found in any species:", paste(missing, collapse = ", ")))
    }
  }
  
  if (length(all_genes) == 0) {
    stop("No genes available.")
  }
  
  # Process each species
  expr_grids <- list()
  for (sp in species) {
    mean_expr <- mean_expr_list[[sp]]
    if (is.null(mean_expr)) {
      warning(paste("No expression data for", sp, ". Setting to 0."))
      full_mean_expr <- matrix(0, nrow = 1, ncol = length(all_genes), dimnames = list("dummy", all_genes))
    } else {
      present_genes <- intersect(all_genes, colnames(mean_expr))
      full_mean_expr <- matrix(0, nrow = nrow(mean_expr), ncol = length(all_genes), dimnames = list(rownames(mean_expr), all_genes))
      full_mean_expr[, present_genes] <- as.matrix(mean_expr[, present_genes])
    }
    expr_grids[[sp]] <- get_scaled_binned_expr(full_mean_expr, sp, num_bins)
  }
  
  # Combine matrices side-by-side
  expr_mat <- do.call(cbind, expr_grids)
  return(expr_mat)
}

# WGCNA Analysis Setup
# Prepare data for WGCNA: transpose expression matrix so rows are "samples" (bins), columns are genes.
#GOI_expr <- Retrive_expr_mat(setdiff(colnames(datExpr), names(net$colors)[net$colors %in% c("brown","green","grey")])) # Use all common genes; optionally pass genes as argument
GOI_expr <- Retrive_expr_mat(GOI)
datExpr <- as.matrix(t(GOI_expr)) # Samples (species_bin) as rows, genes as columns
samples <- rownames(datExpr)
colnames(datExpr)

# Parse sample names to extract species and bin numbers
split_samples <- strsplit(samples, "_")
species_list <- sapply(split_samples, `[`, 1)
bin_list <- as.numeric(sapply(split_samples, `[`, 2))
# Define scaled stages based on the common grid
num_bins <- length(unique(bin_list)) # Should be 30 or whatever was used
grid <- seq(0, 1, length.out = num_bins)
stage_list <- rep(grid, length(unique(species_list))) # Repeat grid for each species
# Create traits data frame
datTraits <- data.frame(
  species = species_list,
  stage = stage_list,
  row.names = samples
)
# Choose soft-thresholding power for WGCNA network construction
powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
power <- sft$powerEstimate
if (is.na(power)) power <- 6 # Default power if no good fit
# Construct co-expression network and identify modules using blockwiseModules
net <- blockwiseModules(
  datExpr,
  power = power,
  TOMType = "signed Nowick 2", # Allows negative TOM for opposites; stricter than "signed Nowick 2"
  minModuleSize = 100, # Smaller for finer separation (with 300 genes)
  mergeCutHeight = 0.15, # Lower to merge less, preserving distinct dynamics
  #networkType = "signed", # Preserves signs for dynamics
  #consensusQuantile = 0, # Minimum TOM across sets (strict consensus)
  verbose = 3
)
# Compute module eigengenes (summaries of module expression profiles)
MEs <- moduleEigengenes(datExpr, net$colors)$eigengenes
colnames(MEs)
plot_me_trends("grey", gene_wise = TRUE, alpha = 0.2, alphacurve = 0.2, highlight_gene = NULL, MEs = MEs, datTraits = datTraits, species_colors = species_colors, net = net, datExpr = datExpr)
plot_me_trends("grey", gene_wise = F, alpha = 0.15, alphacurve = 0.2, highlight_gene = NULL, MEs = MEs, datTraits = datTraits, species_colors = species_colors, net = net, datExpr = datExpr)

plot_dendrogram(net)
plot_me_trends(NULL, colnumber = 2, MEs = MEs, datTraits = datTraits, species_colors = species_colors)

# Create species-specific stage traits for correlation analysis
# Example module queries and operations (adapt as needed)
module <- "brown"
module <- "grey"
module <- "yellow"
module <- "turquoise"
module <- "blue"
module <- "green"
module <- "red"
sort(names(net$colors)[net$colors == module])
length(names(net$colors)[net$colors == module])
mouse_tf <- readLines("mouse TFs.txt")
intersect(names(net$colors)[net$colors == module], mouse_tf)
write.csv(names(net$colors)[net$colors == module], paste0("module_", module, ".csv"))




##############################################################################
# Function to get the module a gene belongs to
# Given a gene name, returns the module color it belongs to (or NA if not found)
get_module_for_gene <- function(gene, net) {
  if (!gene %in% names(net$colors)) {
    warning(paste("Gene", gene, "not found in the network."))
    return(NA)
  }
  module_label <- net$colors[gene]
  return(module_label)
}

# Usage example
get_module_for_gene("Ccnd1", net)

###############################################################################
# Plot module eigengene trends 
# This function plots smoothed trends of module eigengenes over scaled stages for each species.
# If module is NULL (default), it will plot trends for all modules in MEs, arranged in a grid using cowplot.
plot_me_trends <- function(module = NULL, colnumber = 3, gene_wise = FALSE, alpha = 0.05, highlight_gene = NULL, MEs, datTraits, species_colors, net = NULL, datExpr = NULL, alphacurve = 1) {
  unique_species <- unique(datTraits$species)
  if (gene_wise) {
    if (is.null(module)) {
      stop("gene_wise=TRUE requires a specific module to be provided (module cannot be NULL).")
    }
    if (is.null(net) || is.null(datExpr)) {
      stop("gene_wise=TRUE requires net and datExpr to be provided.")
    }
    
    # Make module lowercase for case-insensitivity
    module <- tolower(module)
    available_modules <- gsub("^ME", "", colnames(MEs))
    if (!module %in% available_modules) {
      stop(paste("Module", module, "not found in MEs. Available modules:", paste(sort(available_modules), collapse = ", ")))
    }
    
    # Get module genes
    module_genes <- names(net$colors)[net$colors == module]
    if (length(module_genes) == 0) {
      stop(paste("No genes found in module", module))
    }
    
    # Extract expression data for module genes
    expr_data <- datExpr[, module_genes, drop = FALSE]
    
    # Prepare long format data: stage, species, gene, expression
    plot_data <- data.frame(
      stage = datTraits$stage,
      species = datTraits$species
    )
    plot_data <- cbind(plot_data, as.data.frame(expr_data))
    plot_data_long <- reshape2::melt(plot_data, id.vars = c("stage", "species"), variable.name = "gene", value.name = "expression")
    
    # Reorder species as factors
    plot_data_long$species <- factor(plot_data_long$species, levels = unique_species)
    
    # Determine y-axis limits from data range
    y_limits <- range(plot_data_long$expression, na.rm = TRUE)
    
    # Define color codes
    background_colors <- species_colors
    highlight_color <- "#FAD32F"
    highlight_colors <- setNames(rep(highlight_color, length(unique_species)), unique_species)
    
    # Base plot (no data to allow separate layers)
    p <- ggplot() +
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
      ) +
      coord_cartesian(ylim = c(-3.5,3.5)) +
      labs(
        title = paste("Gene Expression Trends for Module", module),
        subtitle = paste("Across Species (Individual Genes)"),
        x = "Scaled Developmental Stage (0 to 1)",
        y = "Gene Expression",
        color = "Species",
        fill = "Species"
      )
    
    # Add background genes (non-highlighted) with lighter colors and low alpha
    background_data <- plot_data_long[!plot_data_long$gene %in% highlight_gene, ]
    if (nrow(background_data) > 0) {
      for (sp in unique_species) {
        bd_sp <- background_data[background_data$species == sp, ]
        if (nrow(bd_sp) > 0) {
          p <- p + geom_smooth(data = bd_sp, aes(x = stage, y = expression, group = gene),
                               color = alpha(background_colors[sp], alphacurve), fill = alpha(background_colors[sp], alphacurve),
                               method = "loess", se = TRUE, alpha = alpha, size = 0.5, span = 1.5, level = 0.5)
        }
      }
    }
    
    # Add highlighted genes with darker colors, full alpha, thicker lines, and labels
    if (!is.null(highlight_gene)) {
      highlight_data <- plot_data_long[plot_data_long$gene %in% highlight_gene, ]
      if (nrow(highlight_data) == 0) {
        warning("No highlight genes found in the module.")
      } else {
        p <- p + geom_smooth(data = highlight_data, aes(x = stage, y = expression, group = gene, color = species, fill = species),
                             method = "loess", se = TRUE, alpha = 0.5, size = 1.2, span = 1.5, level = 0.5) +
          scale_color_manual(values = highlight_colors) +
          scale_fill_manual(values = highlight_colors)
        
        # Add labels at the end of the curves (stage = 1)
        label_data <- highlight_data[highlight_data$stage == 1, ]
        label_data$label <- paste(label_data$species, label_data$gene)
        p <- p + geom_text(data = label_data, aes(x = stage, y = expression, label = label, color = species),
                           hjust = 1.05, vjust = 0.5, size = 4, fontface = "bold")
      }
    }
    
    return(p)
  }
  
  # Non-gene_wise mode (original behavior)
  if (is.null(module)) {
    modules <- gsub("^ME", "", colnames(MEs))
    plot_list <- lapply(modules, function(mod) {
      # Get the base plot for the module
      me_col <- paste0("ME", mod)
      plot_data <- data.frame(
        ME = MEs[, me_col],
        species = datTraits$species,
        stage = datTraits$stage
      )
      
      # Reorder species as factors for consistent plotting
      plot_data$species <- factor(plot_data$species, levels = unique_species)
      
      # Determine y-axis limits from data range
      y_limits <- range(plot_data$ME, na.rm = TRUE)
      
      p <- ggplot(plot_data, aes(x = stage, y = ME, color = species, group = species)) +
        geom_smooth(aes(fill = species), method = "loess", se = TRUE, alpha = 0.2, size = 0.8, span = 1.5, level = 0.5) +
        scale_color_manual(values = species_colors) +
        scale_fill_manual(values = species_colors) +
        coord_cartesian(ylim = y_limits) +
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
          plot.subtitle = element_blank(),
          axis.title = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom",
          legend.title = element_text(face = "bold"),
          strip.background = element_rect(fill = "grey90", color = "black"),
          strip.text = element_text(face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      
      # Modify for grid: remove x ticks, x title, y title, legend; set title to module name
      p <- p + theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
        coord_cartesian(ylim = c(-0.3,0.3)) +
        labs(title = mod)
      
      return(p)
    })
    
    # Arrange plots in a grid
    arranged_plot <- plot_grid(plotlist = plot_list, ncol = colnumber)
    return(arranged_plot)
  }
  
  # Proceed with single module eigengene plotting
  module <- tolower(module)
  available_modules <- gsub("^ME", "", colnames(MEs))
  if (!module %in% available_modules) {
    stop(paste("Module", module, "not found in MEs. Available modules:", paste(sort(available_modules), collapse = ", ")))
  }
  
  me_col <- paste0("ME", module)
  
  plot_data <- data.frame(
    ME = MEs[, me_col],
    species = datTraits$species,
    stage = datTraits$stage
  )
  
  # Reorder species as factors for consistent plotting
  plot_data$species <- factor(plot_data$species, levels = unique_species)
  
  # Determine y-axis limits from data range
  y_limits <- range(plot_data$ME, na.rm = TRUE)
  
  p <- ggplot(plot_data, aes(x = stage, y = ME, color = species, group = species)) +
    geom_smooth(aes(fill = species), method = "loess", se = TRUE, alpha = 0.05, size = 0.8, span = 1.5, level = 0.5) +
    scale_color_manual(values = species_colors) +
    scale_fill_manual(values = species_colors) +
    coord_cartesian(ylim = y_limits) +
    labs(
      title = paste("Module Eigengene Trends for Module", module),
      subtitle = "Across Species",
      x = "Scaled Developmental Stage (0 to 1)",
      y = "Module Eigengene Expression",
      color = "Species",
      fill = "Species"
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
# Usage examples
species_colors <- c("human" = "#9467BD", "macaque" = "#A87B20", "mouse" = "#2091A8", "ferret" = "#2CA02C", "rat" = "#A83720")
# Plot for a single module
plot_me_trends("turquoise", MEs = MEs, datTraits = datTraits, species_colors = species_colors)
# Plot for a single module with all gene being plotted
plot_me_trends("turquoise", gene_wise = TRUE, alpha = 0.03, alphacurve = 0.01, highlight_gene = NULL, MEs = MEs, datTraits = datTraits, species_colors = species_colors, net = net, datExpr = datExpr)
plot_me_trends("blue", gene_wise = TRUE, alpha = 0.03, alphacurve = 0.03, highlight_gene = NULL, MEs = MEs, datTraits = datTraits, species_colors = species_colors, net = net, datExpr = datExpr)
plot_me_trends("brown", gene_wise = TRUE, alpha = 0.03, alphacurve = 0.3, highlight_gene = NULL, MEs = MEs, datTraits = datTraits, species_colors = species_colors, net = net, datExpr = datExpr)
plot_me_trends("yellow", gene_wise = TRUE, alpha = 0.03, alphacurve = 0.1, highlight_gene = NULL, MEs = MEs, datTraits = datTraits, species_colors = species_colors, net = net, datExpr = datExpr)
plot_me_trends("green", gene_wise = TRUE, alpha = 0.03, alphacurve = 0.1, highlight_gene = NULL, MEs = MEs, datTraits = datTraits, species_colors = species_colors, net = net, datExpr = datExpr)
plot_me_trends("grey", gene_wise = TRUE, alpha = 0.08, alphacurve = 0.2, highlight_gene = NULL, MEs = MEs, datTraits = datTraits, species_colors = species_colors, net = net, datExpr = datExpr)
plot_me_trends("grey", gene_wise = F, alpha = 0.08, alphacurve = 0.2, highlight_gene = NULL, MEs = MEs, datTraits = datTraits, species_colors = species_colors, net = net, datExpr = datExpr)

# Plot for all modules arranged in a grid with 6 columns
plot_me_trends(NULL, colnumber = 2, MEs = MEs, datTraits = datTraits, species_colors = species_colors)


#########################################################################################
# Redesigned plot_coexpression_network function (adapted to new data structure, uses mouse_tf)
# This function visualizes the top connections in a module's co-expression network, highlighting TFs.
plot_coexpression_network <- function(module, datExpr, power, mouse_tf, top_connections = 50, label_top_n = 10) {
  module_genes <- names(net$colors)[net$colors == module]
  TOM <- TOMsimilarityFromExpr(datExpr[, module_genes], power = power)
  
  # Select top connections based on TOM values
  TOM_upper <- TOM[upper.tri(TOM)]
  threshold <- quantile(TOM_upper, 1 - (top_connections / choose(length(module_genes), 2)))
  links <- which(TOM > threshold, arr.ind = TRUE)
  links <- links[links[,1] < links[,2], ]
  links_df <- data.frame(from = module_genes[links[,1]], to = module_genes[links[,2]])
  
  if (nrow(links_df) == 0) {
    stop("No connections above the threshold. Try increasing top_connections.")
  }
  
  # Create igraph network object
  network <- graph_from_data_frame(links_df, directed = FALSE)
  
  # Get genes in the network
  network_genes <- V(network)$name
  
  # Calculate module eigengenes and membership (kME)
  colors <- rep(module, length(module_genes))
  MEs <- moduleEigengenes(datExpr[, module_genes], colors)$eigengenes
  MM <- cor(datExpr[, module_genes], MEs, use = "p")
  
  # Subset membership to network genes
  MM_network <- MM[network_genes, , drop = FALSE]
  
  # Identify TFs in the module and network
  module_tfs <- intersect(module_genes, mouse_tf)
  tfs_in_network <- module_tfs[module_tfs %in% network_genes]
  
  # Set node attributes: size and color based on kME, labels for TFs and top genes
  V(network)$kME <- MM_network[, 1]
  V(network)$size <- 2 + 10 * (V(network)$kME - min(V(network)$kME)) / (max(V(network)$kME) - min(V(network)$kME))
  V(network)$color <- scales::col_numeric("Blues", domain = V(network)$kME)(V(network)$kME)
  V(network)$label <- NA
  
  # Highlight and label TFs
  if (length(tfs_in_network) > 0) {
    tf_indices <- which(V(network)$name %in% tfs_in_network)
    V(network)$color[tf_indices] <- "red"
    V(network)$label[tf_indices] <- V(network)$name[tf_indices]
  }
  
  # Label top non-TF genes by kME
  top_genes <- module_genes[order(MM[, 1], decreasing = TRUE)]
  top_genes_in_network <- top_genes[top_genes %in% network_genes & !(top_genes %in% tfs_in_network)][1:label_top_n]
  top_idx <- which(V(network)$name %in% top_genes_in_network)
  V(network)$label[top_idx] <- V(network)$name[top_idx]
  
  # Create plot using ggraph
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
# For mouse_tf data: Download the list from https://resources.aertslab.org/cistarget/tf_lists/allTFs_mm.txt (one gene symbol per line).
# Alternatively, use AnimalTFDB[](http://bioinfo.life.hust.edu.cn/AnimalTFDB4/#/download) for mouse TF gene symbols.
# Usage example (after downloading and saving the file):
mouse_tf <- readLines("mouse TFs.txt")
plot_coexpression_network("red", datExpr, power, mouse_tf, top_connections = 500, label_top_n = 100)
plot_coexpression_network("turquoise", datExpr, power, mouse_tf, top_connections = 400, label_top_n = 100)
plot_coexpression_network("blue", datExpr, power, mouse_tf, top_connections = 500, label_top_n = 100)
plot_coexpression_network("brown", datExpr, power, mouse_tf, top_connections = 500, label_top_n = 100)
plot_coexpression_network("yellow", datExpr, power, mouse_tf, top_connections = 500, label_top_n = 100)
plot_coexpression_network("green", datExpr, power, mouse_tf, top_connections = 500, label_top_n = 100)
##################################################################
library(dendextend)
plot_dendrogram <- function(net, label_module = NULL) {
  # Get the dendrogram and module assignments (assuming single block for simplicity)
  dend <- net$dendrograms[[1]]
  groups <- net$colors[net$blockGenes[[1]]]
  moduleColors <- groups
  
  # If label_module is provided, modify only the color bar for highlighting
  if (!is.null(label_module)) {
    # Update the color bar to pink for highlighted modules
    moduleColors[moduleColors %in% label_module] <- "black"
  }
  
  # Plot the dendrogram with module colors (branches remain original)
  plotDendroAndColors(
    dend,
    moduleColors,
    "Module Colors",
    dendroLabels = FALSE,
    hang = 0.2,
    addGuide = TRUE,
    guideHang = 0.05,
    main = "Gene Dendrogram and Module Assignments",
    cex.main = 1.2,
    font.main = 2
  )
}
# Usage examples
plot_dendrogram(net) # Standard plot
plot_dendrogram(net, label_module = c("brown"))
###################################################################
###################################################################
# New function to plot expression trends for specified genes
# Plots smoothed curves over scaled stages for each gene in mouse and rat.
# Labels each curve with "species Gene" at the end of the curve.
# Uses species_colors for curve colors.
# Parameters: genes (vector of gene names), datExpr (expression matrix), datTraits (traits with species and stage), species_colors (named vector for mouse and rat).
plot_specified_genes <- function(genes, datExpr, datTraits, species_colors, set_normalize = FALSE) {
  unique_species <- unique(datTraits$species)
  if (length(genes) == 0) {
    stop("No genes provided.")
  }
  
  # Check if all genes are in datExpr
  missing_genes <- setdiff(genes, colnames(datExpr))
  if (length(missing_genes) > 0) {
    warning(paste("Missing genes:", paste(missing_genes, collapse = ", ")))
  }
  genes <- intersect(genes, colnames(datExpr))
  if (length(genes) == 0) {
    stop("None of the provided genes found in datExpr.")
    return()
  }
  
  # Extract expression data for specified genes
  expr_data <- datExpr[, genes, drop = FALSE]
  
  # Normalize per species per gene if requested (z-score)
  if (set_normalize) {
    for (sp in unique_species) {
      sp_rows <- which(datTraits$species == sp)
      for (g in genes) {
        vals <- expr_data[sp_rows, g]
        if (sd(vals, na.rm = TRUE) == 0 || all(is.na(vals))) {
          expr_data[sp_rows, g] <- 0  # Handle constant or all-NA cases
        } else {
          expr_data[sp_rows, g] <- (vals - mean(vals, na.rm = TRUE)) / sd(vals, na.rm = TRUE)
        }
      }
    }
  }
  
  # Prepare long format data: stage, species, gene, expression
  plot_data <- data.frame(
    stage = datTraits$stage,
    species = datTraits$species
  )
  plot_data <- cbind(plot_data, as.data.frame(expr_data))
  plot_data_long <- reshape2::melt(plot_data, id.vars = c("stage", "species"), variable.name = "gene", value.name = "expression")
  
  # Reorder species as factors
  plot_data_long$species <- factor(plot_data_long$species, levels = unique_species)
  
  # Determine y-axis limits from data range
  y_limits <- range(plot_data_long$expression, na.rm = TRUE)
  
  # Set y-label based on normalization
  y_label <- ifelse(set_normalize, "Normalized Gene Expression (Z-score)", "Gene Expression")
  
  # Base plot
  p <- ggplot(plot_data_long, aes(x = stage, y = expression, color = species, fill = species, group = interaction(species, gene))) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.05, size = 1, span = 3, level = 0.5) +
    scale_color_manual(values = species_colors) +
    scale_fill_manual(values = species_colors) +
    coord_cartesian(ylim = y_limits) +
    labs(
      title = "Gene Expression Trends for Specified Genes",
      subtitle = "Across Species",
      x = "Scaled Developmental Stage (0 to 1)",
      y = y_label,
      color = "Species",
      fill = "Species"
    ) +
    theme_minimal(base_size = 14) +
    coord_cartesian(ylim = c(-3,3)) +
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
  
  # Add labels at the end of the curves (stage = 1)
  label_data <- plot_data_long[plot_data_long$stage == max(plot_data_long$stage), ]
  label_data$label <- paste(label_data$species, label_data$gene)
  p <- p + geom_text(data = label_data, aes(label = label), hjust = 0.8, vjust = -1, size = 4, fontface = "bold")
  
  return(p)
}
# Usage example
species_colors <- c("human" = "#9467BD", "macaque" = "#A87B20", "mouse" = "#2091A8", "ferret" = "#2CA02C", "rat" = "#A83720")

plot_specified_genes(c("AXIN2"), datExpr, datTraits, species_colors, set_normalize = TRUE)
plot_specified_genes(c("WNT7B"), datExpr, datTraits, species_colors, set_normalize = T)
plot_specified_genes(c("WNT7B"), datExpr, datTraits, species_colors, set_normalize = T)
plot_specified_genes(c("WLS"), datExpr, datTraits, species_colors, set_normalize = T)
plot_specified_genes(c("LMO2"), datExpr, datTraits, species_colors, set_normalize = T)
plot_specified_genes(c("FRZB"), datExpr, datTraits, species_colors, set_normalize = T)

plot_specified_genes(c("AXIN2", "WNT3","WNT4","WNT5A","WNT5B","WNT7A","WNT7B", "WNT8B", "FRZB", "WLS", "LMO2"), datExpr, datTraits, species_colors)
plot_specified_genes(c("AXIN2"), datExpr, datTraits, species_colors)
plot_specified_genes(c("WLS"), datExpr, datTraits, species_colors)
plot_specified_genes(c("WNT7B"), datExpr, datTraits, species_colors)
plot_specified_genes(c("HMGA2"), datExpr, datTraits, species_colors)
plot_specified_genes(c("CCND1"), datExpr, datTraits, species_colors)
