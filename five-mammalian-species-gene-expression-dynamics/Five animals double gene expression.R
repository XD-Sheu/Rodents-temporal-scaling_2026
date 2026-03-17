library(Seurat)
library(ggplot2)
library(rlang)
library(WGCNA)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(igraph)
library(ggraph)
library(clusterProfiler)
library(org.Hs.eg.db)


########################################################################
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


DefaultAssay(Rat_stemcell_regress) <- "RNA"
Rat_stemcell_regress <- JoinLayers(Rat_stemcell_regress)

DefaultAssay(Mouse_stemcell_regress) <- "RNA"
Mouse_stemcell_regress <- JoinLayers(Mouse_stemcell_regress)

DefaultAssay(Ferret_stemcell_regress) <- "RNA"
Ferret_stemcell_regress <- JoinLayers(Ferret_stemcell_regress)

DefaultAssay(Macaque_stemcell_regress) <- "RNA"
Macaque_stemcell_regress <- JoinLayers(Macaque_stemcell_regress)

DefaultAssay(Human_stemcell_regress) <- "RNA"
Human_stemcell_regress <- JoinLayers(Human_stemcell_regress)


########################################################################
double_gene_expression_plot <- function(object, assay = "RNA", cell_type_metadata, cell_type, 
                                        gene_iden_1, gene_iden_2, x_limit = NULL, y_limit = NULL, 
                                        regress = FALSE, show_axis_titles = TRUE, 
                                        grey_dots = FALSE, regress_nonzero_only = FALSE) {
  require(Seurat)
  require(ggplot2)
  require(rlang)
  
  # Check if the metadata column exists
  if (!cell_type_metadata %in% colnames(object@meta.data)) {
    stop(paste("Metadata column", cell_type_metadata, "not found in object@meta.data"))
  }
  
  # Identify cells belonging to the specified cell type
  cells <- rownames(object@meta.data)[object@meta.data[[cell_type_metadata]] %in% cell_type]
  
  if (length(cells) == 0) {
    stop("No cells found for the specified cell type.")
  }
  
  # Fetch normalized expression data for the two genes
  expr <- FetchData(object, vars = c(gene_iden_1, gene_iden_2), cells = cells, assay = assay, slot = "data")
  
  # Set dot aesthetics based on grey_dots option
  if (grey_dots) {
    dot_color <- "grey"
    dot_alpha <- 0.5
  } else {
    dot_color <- "darkblue"
    dot_alpha <- 0.6
  }
  
  # Create the base scatter plot
  p <- ggplot(expr, aes(x = !!sym(gene_iden_1), y = !!sym(gene_iden_2))) +
    geom_point(color = dot_color, alpha = dot_alpha, size = 1.2) +
    theme_classic() +
    labs(
      x = paste(gene_iden_1, "Expression"),
      y = paste(gene_iden_2, "Expression")
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  # Remove axis titles if requested
  if (!show_axis_titles) {
    p <- p + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  }
  
  # Apply axis limits if provided
  if (!is.null(x_limit) || !is.null(y_limit)) {
    p <- p + coord_cartesian(xlim = x_limit, ylim = y_limit)
  }
  
  # Perform linear regression if requested
  if (regress) {
    # Filter for non-zero cells if requested
    if (regress_nonzero_only) {
      expr_regress <- expr[expr[[gene_iden_1]] > 0 & expr[[gene_iden_2]] > 0, ]
      cat("Regression using", nrow(expr_regress), "cells with non-zero expression in both genes\n")
    } else {
      expr_regress <- expr
    }
    
    # Fit linear model
    model <- lm(as.formula(paste(gene_iden_2, "~", gene_iden_1)), data = expr_regress)
    slope <- coef(model)[2]
    rsq <- summary(model)$r.squared
    
    # Add regression line to plot (using filtered data if applicable)
    p <- p + geom_smooth(data = expr_regress, method = "lm", se = FALSE, color = "red", linewidth = 1)
    
    # Add annotation for slope and R-squared
    annot_text <- paste("Slope: ", round(slope, 3), "\nR²: ", round(rsq, 3), sep = "")
    p <- p + annotate("text", x = Inf, y = Inf, label = annot_text, hjust = 1.1, vjust = 1.1, size = 4, color = "black")
    
    # Output statistics to console
    cat("Linear Regression Results:\n")
    cat("Slope:", round(slope, 3), "\n")
    cat("R-squared:", round(rsq, 3), "\n")
  }
  
  return(p)
}



double_gene_expression_plot(Rat_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "Lhx2", gene_iden_2 = "Lmx1a")
double_gene_expression_plot(Mouse_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "Lhx2", gene_iden_2 = "Lmx1a")
double_gene_expression_plot(Ferret_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "LHX2", gene_iden_2 = "LMX1A")
double_gene_expression_plot(Macaque_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "LHX2", gene_iden_2 = "LMX1A")
double_gene_expression_plot(Human_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "LHX2", gene_iden_2 = "LMX1A")

double_gene_expression_plot(Rat_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "Lhx2", gene_iden_2 = "Emx2")
double_gene_expression_plot(Mouse_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "Lhx2", gene_iden_2 = "Emx2")
double_gene_expression_plot(Ferret_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "LHX2", gene_iden_2 = "EMX2")
double_gene_expression_plot(Macaque_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "LHX2", gene_iden_2 = "EMX2")
double_gene_expression_plot(Human_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "LHX2", gene_iden_2 = "EMX2")

double_gene_expression_plot(Rat_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "Lmx1a", gene_iden_2 = "Wnt2", regress = T, grey_dots = T, regress_nonzero_only = T, show_axis_titles = F)
double_gene_expression_plot(Rat_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "Lmx1a", gene_iden_2 = "Wnt2b", regress = T, grey_dots = T, regress_nonzero_only = T, show_axis_titles = F)
double_gene_expression_plot(Rat_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "Lmx1a", gene_iden_2 = "Wnt3", regress = T, grey_dots = T, regress_nonzero_only = T, show_axis_titles = F)
double_gene_expression_plot(Rat_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "Lmx1a", gene_iden_2 = "Wnt3a", regress = T, grey_dots = T, regress_nonzero_only = T, show_axis_titles = F)
double_gene_expression_plot(Rat_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "Lmx1a", gene_iden_2 = "Wnt4", regress = T, grey_dots = T, regress_nonzero_only = T, show_axis_titles = F)
double_gene_expression_plot(Rat_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "Lmx1a", gene_iden_2 = "Wnt5a", regress = T, grey_dots = T, regress_nonzero_only = T, show_axis_titles = F)
double_gene_expression_plot(Rat_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "Lmx1a", gene_iden_2 = "Wnt5b", regress = T, grey_dots = T, regress_nonzero_only = T, show_axis_titles = F)
double_gene_expression_plot(Rat_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "Lmx1a", gene_iden_2 = "Wnt6", regress = T, grey_dots = T, regress_nonzero_only = T, show_axis_titles = F)
double_gene_expression_plot(Rat_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "Lmx1a", gene_iden_2 = "Wnt7b", regress = T, grey_dots = T, regress_nonzero_only = T, show_axis_titles = F)
double_gene_expression_plot(Rat_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "Lmx1a", gene_iden_2 = "Wnt8b", regress = T, grey_dots = T, regress_nonzero_only = T, show_axis_titles = F)
double_gene_expression_plot(Rat_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "Lmx1a", gene_iden_2 = "Wnt9a", regress = T, grey_dots = T, regress_nonzero_only = T, show_axis_titles = F)
double_gene_expression_plot(Rat_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "Lmx1a", gene_iden_2 = "Wnt16", regress = T, grey_dots = T, regress_nonzero_only = T, show_axis_titles = F)

double_gene_expression_plot(Rat_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "Lhx2", gene_iden_2 = "Wnt2", regress = T, grey_dots = T, regress_nonzero_only = T, show_axis_titles = F)
double_gene_expression_plot(Rat_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "Lhx2", gene_iden_2 = "Wnt2b", regress = T, grey_dots = T, regress_nonzero_only = T, show_axis_titles = F)
double_gene_expression_plot(Rat_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "Lhx2", gene_iden_2 = "Wnt3", regress = T, grey_dots = T, regress_nonzero_only = T, show_axis_titles = F)
double_gene_expression_plot(Rat_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "Lhx2", gene_iden_2 = "Wnt3a", regress = T, grey_dots = T, regress_nonzero_only = T, show_axis_titles = F)
double_gene_expression_plot(Rat_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "Lhx2", gene_iden_2 = "Wnt4", regress = T, grey_dots = T, regress_nonzero_only = T, show_axis_titles = F)
double_gene_expression_plot(Rat_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "Lhx2", gene_iden_2 = "Wnt5a", regress = T, grey_dots = T, regress_nonzero_only = T, show_axis_titles = F)
double_gene_expression_plot(Rat_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "Lhx2", gene_iden_2 = "Wnt5b", regress = T, grey_dots = T, regress_nonzero_only = T, show_axis_titles = F)
double_gene_expression_plot(Rat_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "Lhx2", gene_iden_2 = "Wnt6", regress = T, grey_dots = T, regress_nonzero_only = T, show_axis_titles = F)
double_gene_expression_plot(Rat_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "Lhx2", gene_iden_2 = "Wnt7b", regress = T, grey_dots = T, regress_nonzero_only = T, show_axis_titles = F)
double_gene_expression_plot(Rat_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "Lhx2", gene_iden_2 = "Wnt8b", regress = T, grey_dots = T, regress_nonzero_only = T, show_axis_titles = F)
double_gene_expression_plot(Rat_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "Lhx2", gene_iden_2 = "Wnt9a", regress = T, grey_dots = T, regress_nonzero_only = T, show_axis_titles = F)
double_gene_expression_plot(Rat_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "Lhx2", gene_iden_2 = "Wnt16", regress = T, grey_dots = T, regress_nonzero_only = T, show_axis_titles = F)

########################################################################
ratio_positive_plot <- function(object, assay = "RNA", cell_type_metadata, cell_type, gene_iden_1, gene_iden_2) {
  require(Seurat)
  require(ggplot2)
  require(rlang)
  require(reshape2)
  
  # Check if the metadata column exists
  if (!cell_type_metadata %in% colnames(object@meta.data)) {
    stop(paste("Metadata column", cell_type_metadata, "not found in object@meta.data"))
  }
  
  # Identify cells belonging to the specified cell type
  cells <- rownames(object@meta.data)[object@meta.data[[cell_type_metadata]] %in% cell_type]
  
  if (length(cells) == 0) {
    stop("No cells found for the specified cell type.")
  }
  
  # Fetch normalized expression data for the two genes
  expr <- FetchData(object, vars = c(gene_iden_1, gene_iden_2), cells = cells, assay = assay, slot = "data")
  
  # Calculate positives
  gene2_pos <- expr[[gene_iden_2]] > 0.1
  double_pos <- expr[[gene_iden_1]] > 0.1 & gene2_pos
  num_double <- sum(double_pos)
  num_gene2 <- sum(gene2_pos)
  
  if (num_gene2 == 0) {
    stop("No cells with positive expression for gene_iden_2.")
  }
  
  # Calculate percentages
  double_ratio <- (num_double / num_gene2) * 100
  single_ratio <- 100 - double_ratio
  
  # Create data frame for plotting
  plot_data <- data.frame(
    category = paste(gene_iden_1, "+ among", gene_iden_2, "+"),
    Double_Positive = double_ratio,
    Single_Positive = single_ratio
  )
  
  # Reshape data for stacked bar plot
  plot_data_long <- reshape2::melt(plot_data, id.vars = "category", variable.name = "Expression", value.name = "Percentage")
  
  # Reverse the order of the fill levels to put Double_Positive at the bottom
  plot_data_long$Expression <- factor(plot_data_long$Expression, levels = c("Single_Positive", "Double_Positive"))
  
  # Create the stacked bar plot
  p <- ggplot(plot_data_long, aes(x = category, y = Percentage, fill = Expression)) +
    geom_bar(stat = "identity", width = 0.6) +
    scale_fill_manual(values = c("Double_Positive" = "darkblue", "Single_Positive" = "lightgrey"),
                      labels = c("Double Positive", paste(gene_iden_2, "+ Only"))) +
    theme_classic() +
    labs(
      x = "",
      y = "Percentage (%)",
      title = paste("Ratio of Double Positive Cells (", gene_iden_1, "+ in", gene_iden_2, "+ )"),
      subtitle = paste("Cell Type:", cell_type, "| Assay:", assay, "| Total Cells:", length(cells))
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 0, vjust = 0.5),
      legend.title = element_blank(),
      legend.position = "top"
    ) +
    ylim(0, 100) +
    annotate(
      "text",
      x = 1,
      y = double_ratio + 5,
      label = paste(num_double, "/", num_gene2),
      size = 4,
      color = "black",
      vjust = -0.5
    )
  
  return(p)
}

ratio_positive_plot(Rat_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "Lmx1a", gene_iden_2 = "Lhx2")
ratio_positive_plot(Mouse_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "Lmx1a", gene_iden_2 = "Lhx2")
ratio_positive_plot(Ferret_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "LMX1A", gene_iden_2 = "LHX2")
ratio_positive_plot(Macaque_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "LMX1A", gene_iden_2 = "LHX2")
ratio_positive_plot(Human_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                            gene_iden_1 = "LMX1A", gene_iden_2 = "LHX2")

ratio_positive_plot(Rat_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                    gene_iden_1 = "Emx2", gene_iden_2 = "Lhx2")
ratio_positive_plot(Mouse_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                    gene_iden_1 = "Emx2", gene_iden_2 = "Lhx2")
ratio_positive_plot(Ferret_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                    gene_iden_1 = "EMX2", gene_iden_2 = "LHX2")
ratio_positive_plot(Macaque_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                    gene_iden_1 = "EMX2", gene_iden_2 = "LHX2")
ratio_positive_plot(Human_stemcell_regress, assay = "SCT", cell_type_metadata = "stemcell_type", cell_type = "vRG",
                    gene_iden_1 = "EMX2", gene_iden_2 = "LHX2")
