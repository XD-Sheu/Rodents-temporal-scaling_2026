library(ggplot2)
library(monocle3)
library(dplyr)
library(Seurat)
library(patchwork)
library(SeuratWrappers)
library(scRNAseq)
library(viridis)
library(dittoSeq)
library(SingleR)
library(ggplot2)
library(RColorBrewer)
library(tidyverse )
library(ggpubr)
library(UCell)
library(glmGamPoi)
library(colorRamp2)
library(ComplexHeatmap)
library(DESeq2)
library(circlize)
library(magick)
library(gridExtra)
library(data.table)
library(tidyverse)
library(ggpubr)
library(UCell)
library(glmGamPoi)
library(colorRamp2)
library(ComplexHeatmap)
library(DESeq2)
library(clusterProfiler)
library(wordcloud)
library(enrichplot)
library(ggupset)
library(scplotter)
library(EnhancedVolcano)
library(psupertime)
library(topGO)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(cowplot)
library(data.table)
library(scales)
library(xfun)
library(SCP)
library(reshape2)
library(stats)
library(scales)

setwd("F:/Single cell analysis/Mouse and rat/Network and trend/objects")
data_integrated <- read_rds("data_integrated_25_07_28.rds")

Mouse_stemcell <- subset(subset(data_integrated, subset = species %in% "mouse"), subset = Celltype %in% c("Neuroepithelium cell (NEC)","Radial glia cell (RGC)"))
Rat_stemcell <- subset(subset(data_integrated, subset = species %in% "rat"), subset = Celltype %in% c("Neuroepithelium cell (NEC)","Radial glia cell (RGC)"))


dayl <- c("1","2","3")
names(dayl) <- c("E12", "E14", "E16")
Idents(Mouse_stemcell) <- "day"
Mouse_stemcell <- RenameIdents(Mouse_stemcell, dayl)
Mouse_stemcell$stage_order <- Idents(Mouse_stemcell)

dayl <- c("1","2","3")
names(dayl) <- c("E14", "E16", "E18")
Idents(Rat_stemcell) <- "day"
Rat_stemcell <- RenameIdents(Rat_stemcell, dayl)
Rat_stemcell$stage_order <- Idents(Rat_stemcell)

saveRDS(Mouse_stemcell, "Mouse_stemcell.rds")
saveRDS(Rat_stemcell, "Rat_stemcell.rds")

#Mouse_stemcell <- read_rds("Mouse_stemcell.rds")
#Rat_stemcell <- read_rds("Rat_stemcell.rds")

##########################################################
# Mouse
#Psupertime based trajecotry analysis

# Set metacells to make sure each stage have ~50 metacells
Mouse_stemcell$metacell <- NA

stages <- unique(Mouse_stemcell$stage_order)

for (stage in stages) {
  message("Processing stage: ", stage)
  seurat_subset <- subset(Mouse_stemcell, subset = stage_order == stage)
  seurat_subset <- FindClusters(seurat_subset, resolution = 2, verbose = FALSE)
  
  initial_clusters <- seurat_subset$seurat_clusters
  metacells <- split(colnames(seurat_subset), initial_clusters)
  message("Initial FindClusters produced ", length(metacells), " clusters")
  
  while (length(metacells) < 100) {
    sizes <- sapply(metacells, length)
    if (max(sizes) < 30) {
      message("Cluster with fewer than 30 cells detected. Stopping subclustering for stage ", stage)
      break
    }
    largest_idx <- which.max(sizes)
    largest_metacell <- metacells[[largest_idx]]
    seurat_sub <- subset(seurat_subset, cells = largest_metacell)
    seurat_sub <- FindClusters(seurat_sub, resolution = 1, verbose = FALSE)
    new_clusters <- seurat_sub$seurat_clusters
    new_metacells <- split(largest_metacell, new_clusters)
    message("Subclustering added ", length(new_metacells), " new clusters")
    metacells <- metacells[-largest_idx]
    metacells <- c(metacells, new_metacells)
  }
  metacell_labels <- paste0(stage, "_", 1:length(metacells))
  for (i in 1:length(metacells)) {
    Mouse_stemcell$metacell[metacells[[i]]] <- metacell_labels[i]
  }
  message("Stage ", stage, " completed with ", length(metacells), " metacells")
}

DefaultAssay(Mouse_stemcell) <- "RNA"
averaged_exp <- AggregateExpression(Mouse_stemcell, group.by = "metacell", assays = "RNA", slot = "data")
averaged_exp <- averaged_exp[[1]]
colnames(averaged_exp) <- sub("^g(.*)", "\\1", colnames(averaged_exp))
colnames(averaged_exp) <- gsub("-", "_", colnames(averaged_exp))

metacell_types <- sapply(unique(Mouse_stemcell$metacell), function(m) {
  cells <- which(Mouse_stemcell$metacell == m)
  cell_types <- Mouse_stemcell$Celltype[cells]
  tbl <- table(cell_types)
  names(tbl)[which.max(tbl)]
})

stages <- sapply(strsplit(colnames(averaged_exp), "_"), function(x) x[1])
meta_data <- data.frame(
  stage_order = stages,
  Celltype = metacell_types[colnames(averaged_exp)],
  row.names = colnames(averaged_exp)
)

Mouse_stemcell_meta <- CreateSeuratObject(counts = averaged_exp, meta.data = meta_data, project = "Mouse_stemcell_meta")
Mouse_stemcell_meta <- NormalizeData(Mouse_stemcell_meta, normalization.method = "LogNormalize", scale.factor = 10000)
DefaultAssay(Mouse_stemcell_meta) <- "RNA"
Idents(Mouse_stemcell_meta) <- "stage_order"
stage_order_sorted <- c("1","2","3")
Idents(Mouse_stemcell_meta) <- factor(Idents(Mouse_stemcell_meta), levels= stage_order_sorted)
Mouse_stemcell_meta$stage_order_sort <- Idents(Mouse_stemcell_meta)
Idents(Mouse_stemcell_meta) <- "stage_order_sort"
Mouse_stemcell_meta <- FindVariableFeatures(Mouse_stemcell_meta, selection.method = "vst", nfeatures = 4000)
mouse_hvg <- unique(c(VariableFeatures(Mouse_stemcell),VariableFeatures(Mouse_stemcell_meta)))
setwd("F:/Single cell analysis/Mouse and rat/Network and trend/objects")
saveRDS(Mouse_stemcell_meta,"Mouse_stemcell_meta_25_07_28.rds")
#Mouse_stemcell_meta <- read_rds("Mouse_stemcell_meta_25_07_28.rds")


#OLR model
Mouse_OLR <- as.SingleCellExperiment(Mouse_stemcell_meta, assay = "RNA", slot = "data")
conserved_gene_list <- intersect(rownames(Mouse_stemcell_meta@assays$RNA$data), rownames(Rat_stemcell_meta@assays$RNA$data))
conserved_gene_list <- conserved_gene_list[!grepl("^Gm|LOC|Rik$", conserved_gene_list)]
scoremodel_Mouse <- psupertime(Mouse_OLR, assay_type = "logcounts", Mouse_OLR$stage_order_sort, sel_genes = "list", gene_list = conserved_gene_list, n_folds = 3, method = "proportional")
#scoremodel_Mouse <- psupertime(Mouse_OLR, Mouse_OLR$stage_order_sort, n_folds = 5, method = "proportional")

plot_labels_over_psupertime_homemade(scoremodel_Mouse, palette = "YLGnBu")
plot_identified_genes_over_psupertime_homemade (scoremodel_Mouse, n_to_plot=63)
plot_specified_genes_over_psupertime_homemade(scoremodel_Mouse, Interest_top, palette = "Blues")
plot_specified_genes_over_psupertime_homemade(scoremodel_Mouse, c("Nsd1","Meis2","Auts2", "Akt3", "Ezh2", "Eed", "Hmga2", "Ccnd1"), palette = "Blues")

setwd("F:/Single cell analysis/Mouse and rat/Network and trend/objects")
saveRDS(scoremodel_Mouse, "scoremodel_Mouse_25_07_28.rds")
scoremodel_Mouse <- read_rds("scoremodel_Mouse_25_07_28.rds")



# Perform linear regression for each gene
expr_matrix <- t(scoremodel_Mouse$x_data)
stages <- scoremodel_Mouse$proj_dt
stages <- stages[order(stages$psuper),]
expr_matrix <- expr_matrix[ , stages$cell_id]
stages <- stages$psuper
stage_levels <- unique(stages)

mean_expr <- do.call(rbind, lapply(stage_levels, function(s) {
  cells <- which(stages == s)
  if (length(cells) > 0) {
    rowMeans(expr_matrix[, cells, drop = FALSE])
  } else {
    rep(NA, nrow(expr_matrix))
  }
}))
rownames(mean_expr) <- stage_levels
colnames(mean_expr) <- rownames(expr_matrix)
saveRDS(mean_expr, "mouse_mean_expr.rds")

# Initialize results storage
message("Fitting linear regression for ", ncol(mean_expr), " genes")
results <- data.frame(
  gene = colnames(mean_expr),
  coefficient = NA,
  p_value = NA,
  adj_p_value = NA,
  adj_r_squared = NA,
  stringsAsFactors = FALSE
)

# Fit linear regression for each gene
#stage_numeric <- seq_along(stage_levels)
stage_numeric <- rescale(as.numeric(stage_levels), to = c(0, 1))

for (i in 1:ncol(mean_expr)) {
  expr <- mean_expr[, i]
  if (any(is.na(expr))) next  # Skip genes with missing data
  model <- lm(expr ~ stage_numeric)
  summary_model <- summary(model)
  results$coefficient[i] <- coef(model)[2]  # Slope for stage_numeric
  results$p_value[i] <- summary_model$coefficients[2, 4]  # P-value for stage_numeric
  results$adj_r_squared[i] <- summary_model$adj.r.squared
}

results$adj_p_value <- p.adjust(results$p_value, method = "BH")

# Filter for significant linear trends and rank by adjusted R-squared
results <- results[!is.na(results$adj_r_squared), ]
results <- results[results$adj_p_value < 0.05, ]
top_500_genes <- results[order(-results$adj_r_squared), ][1:500, ]
setwd("F:/Single cell analysis/Mouse and rat/Network and trend/objects")
write.csv(top_500_genes, "Top 500 linear genes Mouse.csv")
write.csv(results, "Linear genes all Mouse.csv")
top_500_genes <- read.csv("Top 500 linear genes Mouse.csv")
top_15_increasing <- top_500_genes[top_500_genes$coefficient > 0, ][1:15, ]$gene
top_15_decreasing <- top_500_genes[top_500_genes$coefficient < 0, ][1:15, ]$gene
Interest_30 <- c(top_15_increasing, top_15_decreasing)
Interest_30
plot_specified_genes_over_psupertime_homemade(scoremodel_Mouse, Interest_30, palette = "Blues")



##########################################################
# Rat

# Set metacells to make sure each stage have ~50 metacells
Rat_stemcell$metacell <- NA

stages <- unique(Rat_stemcell$stage_order)

for (stage in stages) {
  message("Processing stage: ", stage)
  seurat_subset <- subset(Rat_stemcell, subset = stage_order == stage)
  seurat_subset <- FindClusters(seurat_subset, resolution = 2, verbose = FALSE)
  
  initial_clusters <- seurat_subset$seurat_clusters
  metacells <- split(colnames(seurat_subset), initial_clusters)
  message("Initial FindClusters produced ", length(metacells), " clusters")
  
  while (length(metacells) < 100) {
    sizes <- sapply(metacells, length)
    if (max(sizes) < 30) {
      message("Cluster with fewer than 30 cells detected. Stopping subclustering for stage ", stage)
      break
    }
    largest_idx <- which.max(sizes)
    largest_metacell <- metacells[[largest_idx]]
    seurat_sub <- subset(seurat_subset, cells = largest_metacell)
    seurat_sub <- FindClusters(seurat_sub, resolution = 1, verbose = FALSE)
    new_clusters <- seurat_sub$seurat_clusters
    new_metacells <- split(largest_metacell, new_clusters)
    message("Subclustering added ", length(new_metacells), " new clusters")
    metacells <- metacells[-largest_idx]
    metacells <- c(metacells, new_metacells)
  }
  metacell_labels <- paste0(stage, "_", 1:length(metacells))
  for (i in 1:length(metacells)) {
    Rat_stemcell$metacell[metacells[[i]]] <- metacell_labels[i]
  }
  message("Stage ", stage, " completed with ", length(metacells), " metacells")
}

DefaultAssay(Rat_stemcell) <- "RNA"
averaged_exp <- AggregateExpression(Rat_stemcell, group.by = "metacell", assays = "RNA", slot = "data")
averaged_exp <- averaged_exp[[1]]
colnames(averaged_exp) <- sub("^g(.*)", "\\1", colnames(averaged_exp))
colnames(averaged_exp) <- gsub("-", "_", colnames(averaged_exp))

metacell_types <- sapply(unique(Rat_stemcell$metacell), function(m) {
  cells <- which(Rat_stemcell$metacell == m)
  cell_types <- Rat_stemcell$Celltype[cells]
  tbl <- table(cell_types)
  names(tbl)[which.max(tbl)]
})

stages <- sapply(strsplit(colnames(averaged_exp), "_"), function(x) x[1])
meta_data <- data.frame(
  stage_order = stages,
  Celltype = metacell_types[colnames(averaged_exp)],
  row.names = colnames(averaged_exp)
)

Rat_stemcell_meta <- CreateSeuratObject(counts = averaged_exp, meta.data = meta_data, project = "Rat_stemcell_meta")
Rat_stemcell_meta <- NormalizeData(Rat_stemcell_meta, normalization.method = "LogNormalize", scale.factor = 10000)
DefaultAssay(Rat_stemcell_meta) <- "RNA"
Idents(Rat_stemcell_meta) <- "stage_order"
stage_order_sorted <- c("1","2","3")
Idents(Rat_stemcell_meta) <- factor(Idents(Rat_stemcell_meta), levels= stage_order_sorted)
Rat_stemcell_meta$stage_order_sort <- Idents(Rat_stemcell_meta)
Idents(Rat_stemcell_meta) <- "stage_order_sort"
Rat_stemcell_meta <- FindVariableFeatures(Rat_stemcell_meta, selection.method = "vst", nfeatures = 4000)
rat_hvg <- unique(c(VariableFeatures(Rat_stemcell),VariableFeatures(Rat_stemcell_meta)))
setwd("F:/Single cell analysis/Mouse and rat/Network and trend/objects")
saveRDS(Rat_stemcell_meta,"Rat_stemcell_meta_25_07_28.rds")
#Rat_stemcell_meta <- read_rds("Rat_stemcell_meta_25_07_28.rds")
#rodent_hvg <- intersect(mouse_hvg, rat_hvg)
#saveRDS(Rodent_hvg, "rodent_hvg_gene.rds")


#OLR model
Rat_OLR <- as.SingleCellExperiment(Rat_stemcell_meta, assay = "RNA")
conserved_gene_list <- intersect(rownames(Mouse_stemcell_meta@assays$RNA$data), rownames(Rat_stemcell_meta@assays$RNA$data))
conserved_gene_list <- conserved_gene_list[!grepl("^Gm|Rik$", conserved_gene_list)]
scoremodel_Rat <- psupertime(Rat_OLR, assay_type = "logcounts", Rat_OLR$stage_order_sort, sel_genes = "list", gene_list = conserved_gene_list, n_folds = 3, method = "proportional")
#scoremodel_Rat <- psupertime(Rat_OLR, Rat_OLR$stage_order_sort, n_folds = 5, method = "proportional")

plot_labels_over_psupertime_homemade(scoremodel_Rat, palette = "YLGnBu")
plot_identified_genes_over_psupertime_homemade (scoremodel_Rat, n_to_plot=63)
plot_specified_genes_over_psupertime_homemade(scoremodel_Rat, Interest_top, palette = "Blues")
plot_specified_genes_over_psupertime_homemade(scoremodel_Rat, c("NSD1","MEIS2","AUTS2", "AKT3", "EZH2", "EED", "HMGA2"), palette = "Blues")
view(scoremodel_Rat$beta_dt)
view(goanalysis$clusters_dt)
view(scoremodel$x_data)

saveRDS(scoremodel_Rat, "scoremodel_Rat_25_07_28.rds")
scoremodel_Rat <- read_rds("scoremodel_Rat_25_07_28.rds")



# Perform linear regression for each gene
expr_matrix <- t(scoremodel_Rat$x_data)
stages <- scoremodel_Rat$proj_dt
stages <- stages[order(stages$psuper),]
expr_matrix <- expr_matrix[ , stages$cell_id]
stages <- stages$psuper
stage_levels <- unique(stages)

mean_expr <- do.call(rbind, lapply(stage_levels, function(s) {
  cells <- which(stages == s)
  if (length(cells) > 0) {
    rowMeans(expr_matrix[, cells, drop = FALSE])
  } else {
    rep(NA, nrow(expr_matrix))
  }
}))
rownames(mean_expr) <- stage_levels
colnames(mean_expr) <- rownames(expr_matrix)
saveRDS(mean_expr, "rat_mean_expr.rds")

# Initialize results storage
message("Fitting linear regression for ", ncol(mean_expr), " genes")
results <- data.frame(
  gene = colnames(mean_expr),
  coefficient = NA,
  p_value = NA,
  adj_p_value = NA,
  adj_r_squared = NA,
  stringsAsFactors = FALSE
)

# Fit linear regression for each gene
#stage_numeric <- seq_along(stage_levels)
stage_numeric <- rescale(as.numeric(stage_levels), to = c(0, 1))

for (i in 1:ncol(mean_expr)) {
  expr <- mean_expr[, i]
  if (any(is.na(expr))) next  # Skip genes with missing data
  model <- lm(expr ~ stage_numeric)
  summary_model <- summary(model)
  results$coefficient[i] <- coef(model)[2]  # Slope for stage_numeric
  results$p_value[i] <- summary_model$coefficients[2, 4]  # P-value for stage_numeric
  results$adj_r_squared[i] <- summary_model$adj.r.squared
}

results$adj_p_value <- p.adjust(results$p_value, method = "BH")

# Filter for significant linear trends and rank by adjusted R-squared
results <- results[!is.na(results$adj_r_squared), ]
results <- results[results$adj_p_value < 0.05, ]
top_500_genes <- results[order(-results$adj_r_squared), ][1:500, ]
write.csv(top_500_genes, "Top 500 linear genes Rat.csv")
write.csv(results, "Linear genes all Rat.csv")
top_15_increasing <- top_500_genes[top_500_genes$coefficient > 0, ][1:15, ]$gene
top_15_decreasing <- top_500_genes[top_500_genes$coefficient < 0, ][1:15, ]$gene
Interest_30 <- c(top_15_increasing, top_15_decreasing)
Interest_30
plot_specified_genes_over_psupertime_homemade(scoremodel_Rat, Interest_30, palette = "Blues")



