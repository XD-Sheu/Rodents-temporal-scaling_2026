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

setwd("F:/Single cell analysis/Parallel development/Human/Objects")
Human_stemcell_regress <- read_rds("Human_stemcell_regress_25_04_29.rds")
Human_vRG <- subset(subset(Human_stemcell_regress, subset = Estimated_postconceptional_age_in_days %in% c("54","60","63","77","84")), subset = stemcell_type %in% c("vRG"))
Human_ExN_IPC <- subset(subset(Human_stemcell_regress, subset = Estimated_postconceptional_age_in_days %in% c("54","60","63","77","84")), subset = stemcell_type %in% c("ExN-IPC"))
Human_oRG <- subset(subset(Human_stemcell_regress, subset = Estimated_postconceptional_age_in_days %in% c("54","60","63","77","84")), subset = stemcell_type %in% c("oRG"))

setwd("F:/Single cell analysis/Parallel development/Macaque/Objects")
Macaque_stemcell_regress <- read_rds("Macaque_stemcell_regress_25_05_19.rds")
Macaque_vRG <- subset(subset(Macaque_stemcell_regress, subset = age %in% c("E42","E54","E64")), subset = stemcell_type %in% c("vRG"))
Macaque_ExN_IPC <- subset(subset(Macaque_stemcell_regress, subset = age %in% c("E42","E54","E64")), subset = stemcell_type %in% c("ExN-IPC"))
Macaque_oRG <- subset(subset(Macaque_stemcell_regress, subset = age %in% c("E42","E54","E64")), subset = stemcell_type %in% c("oRG"))

setwd("F:/Single cell analysis/Parallel development/Ferret/Objects")
Ferret_stemcell_regress <- read_rds("Ferret_stemcell_regress_25_05_19.rds")
Ferret_vRG <- subset(subset(Ferret_stemcell_regress, subset = stage %in% c("E25","E34","E40")), subset = stemcell_type %in% c("vRG"))
Ferret_ExN_IPC <- subset(subset(Ferret_stemcell_regress, subset = stage %in% c("E25","E34","E40")), subset = stemcell_type %in% c("ExN-IPC"))
Ferret_oRG <- subset(subset(Ferret_stemcell_regress, subset = stage %in% c("E25","E34","E40")), subset = stemcell_type %in% c("oRG"))

setwd("F:/Single cell analysis/Parallel development/Mouse/Objects")
Mouse_stemcell_regress <- read_rds("Mouse_stemcell_regress_25_05_19.rds")
Mouse_vRG <- subset(subset(Mouse_stemcell_regress, subset = stage %in% c("E12","E13","E14","E15","E16")), subset = stemcell_type %in% c("vRG"))
Mouse_ExN_IPC <- subset(subset(Mouse_stemcell_regress, subset = stage %in% c("E12","E13","E14","E15","E16")), subset = stemcell_type %in% c("ExN-IPC"))
Mouse_oRG <- subset(subset(Mouse_stemcell_regress, subset = stage %in% c("E12","E13","E14","E15","E16")), subset = stemcell_type %in% c("oRG"))

setwd("F:/Single cell analysis/Mouse and rat/Network and trend/objects")
Rat_stemcell_regress <- read_rds("Rat_stemcell_regress_25_08_06.rds")
dayl <- c("1","2","3")
names(dayl) <- c("E14", "E16", "E18")
Idents(Rat_stemcell_regress) <- "day"
Rat_stemcell_regress <- RenameIdents(Rat_stemcell_regress, dayl)
Rat_stemcell_regress$stage_order <- Idents(Rat_stemcell_regress)
Rat_vRG <- subset(subset(Rat_stemcell_regress, subset = stage_order %in% c("1","2","3")), subset = stemcell_type %in% c("vRG", "oRG"))
Rat_ExN_IPC <- subset(subset(Rat_stemcell_regress, subset = stage_order %in% c("1","2","3")), subset = stemcell_type %in% c("ExN-IPC"))
Rat_oRG <- subset(subset(Rat_stemcell_regress, subset = stage_order %in% c("1","2","3")), subset = stemcell_type %in% c("oRG"))

Gene_atlas <- read.csv("Atlas conserved age genes with linear model and HAR and disease.csv", header = T)
Commongene <- intersect(Gene_atlas$Conserved.genes, toupper(rownames(Rat_stemcell_regress)))
saveRDS(Commongene, "Five_animal_shared_genes.rds")
Commongene <- read_rds("Five_animal_shared_genes.rds")


saveRDS(Human_vRG, "Human_vRG.rds")
saveRDS(Human_ExN_IPC, "Human_ExN_IPC.rds")
saveRDS(Human_oRG, "Human_oRG.rds")
saveRDS(Macaque_vRG, "Macaque_vRG.rds")
saveRDS(Macaque_ExN_IPC, "Macaque_ExN_IPC.rds")
saveRDS(Macaque_oRG, "Macaque_oRG.rds")
saveRDS(Ferret_vRG, "Ferret_vRG.rds")
saveRDS(Ferret_ExN_IPC, "Ferret_ExN_IPC.rds")
saveRDS(Ferret_oRG, "Ferret_oRG.rds")
saveRDS(Mouse_vRG, "Mouse_vRG.rds")
saveRDS(Mouse_ExN_IPC, "Mouse_ExN_IPC.rds")
saveRDS(Mouse_oRG, "Mouse_oRG.rds")
saveRDS(Rat_vRG, "Rat_vRG.rds")
saveRDS(Rat_ExN_IPC, "Rat_ExN_IPC.rds")
saveRDS(Rat_oRG, "Rat_oRG.rds")
#Human_vRG <- read_rds("Human_vRG.rds")
#Macaque_vRG <- read_rds("Macaque_vRG.rds")
#Ferret_vRG <- read_rds("Ferret_vRG.rds")
#Mouse_vRG <- read_rds("Mouse_vRG.rds")
#Rat_vRG <- read_rds("Rat_vRG.rds")

##########################################################
get_mean_expr_after_OLR <- function(seurat_obj, metacell_num = 100, set_stage_order) {
  seurat_obj$metacell <- NA
  Idents(seurat_obj) <- "stage_order"
  unique_stages <- unique(Idents(seurat_obj))
  sorted_stages_num <- sort(as.numeric(as.character(unique_stages)))
  remapped_levels <- as.character(seq_along(sorted_stages_num))
  sorted_stages_char <- as.character(sorted_stages_num)
  mapping <- setNames(as.character(remapped_levels), sorted_stages_char)
  seurat_obj <- RenameIdents(seurat_obj, mapping)
  seurat_obj$stage_order <- Idents(seurat_obj)
  Idents(seurat_obj) <- "stage_order"
  
  stages <- unique(seurat_obj$stage_order)
  for (stage in stages) {
    message("Processing stage: ", stage)
    seurat_subset <- subset(seurat_obj, subset = stage_order %in% as.character(stage))
    seurat_subset <- FindClusters(seurat_subset, resolution = 2, verbose = FALSE)
    
    initial_clusters <- seurat_subset$seurat_clusters
    metacells <- split(colnames(seurat_subset), initial_clusters)
    message("Initial FindClusters produced ", length(metacells), " clusters")
    
    while (length(metacells) < metacell_num) {
      sizes <- sapply(metacells, length)
      if (max(sizes) < 50) {
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
      seurat_obj$metacell[metacells[[i]]] <- metacell_labels[i]
    }
    message("Stage ", stage, " completed with ", length(metacells), " metacells")
  }
  DefaultAssay(seurat_obj) <- "SCT"
  averaged_exp <- AggregateExpression(seurat_obj, group.by = "metacell", assays = "SCT", slot = "data")
  averaged_exp <- averaged_exp[[1]]
  colnames(averaged_exp) <- sub("^g(.*)", "\\1", colnames(averaged_exp))
  colnames(averaged_exp) <- gsub("-", "_", colnames(averaged_exp))
  rownames(averaged_exp) <- sub("\\.BLAST.*", "", rownames(averaged_exp))
  rownames(averaged_exp) <- toupper(rownames(averaged_exp))
  metacell_types <- sapply(unique(seurat_obj$metacell), function(m) {
    cells <- which(seurat_obj$metacell == m)
    cell_types <- seurat_obj$stemcell_type[cells]
    tbl <- table(cell_types)
    names(tbl)[which.max(tbl)]
  })
  stages <- sapply(strsplit(colnames(averaged_exp), "_"), function(x) x[1])
  meta_data <- data.frame(
    stage_order = stages,
    Celltype = metacell_types[colnames(averaged_exp)],
    row.names = colnames(averaged_exp)
  )
  rownames(averaged_exp) <- make.unique(rownames(averaged_exp))
  Stemcell_meta <- CreateSeuratObject(counts = averaged_exp, meta.data = meta_data, project = "Stemcell_meta")
  Stemcell_meta <- NormalizeData(Stemcell_meta, normalization.method = "LogNormalize", scale.factor = 10000)
  Stemcell_meta <- FindVariableFeatures(Stemcell_meta, selection.method = "vst", nfeatures = 5000)
  Hvg <- VariableFeatures(Stemcell_meta)
  #saveRDS(Hvg, "rat_vRG_HVG.rds")
  DefaultAssay(Stemcell_meta) <- "RNA"
  Idents(Stemcell_meta) <- "stage_order"
  
  # Dynamically remap stage_order to consecutive integers starting from 1
  unique_stages <- unique(Idents(Stemcell_meta))
  sorted_stages_num <- sort(as.numeric(as.character(unique_stages)))
  remapped_levels <- as.character(seq_along(sorted_stages_num))
  sorted_stages_char <- as.character(sorted_stages_num)
  mapping <- setNames(as.character(remapped_levels), sorted_stages_char)
  Stemcell_meta <- RenameIdents(Stemcell_meta, mapping)
  Stemcell_meta$stage_order_sort <- Idents(Stemcell_meta)
  Idents(Stemcell_meta) <- "stage_order_sort"

  stage_order_sorted <- set_stage_order
  Idents(Stemcell_meta) <- factor(Idents(Stemcell_meta), levels= set_stage_order)
  Stemcell_meta$stage_order_sort <- Idents(Stemcell_meta)
  Idents(Stemcell_meta) <- "stage_order_sort"
  
  
  # OLR model
  Species_OLR <- as.SingleCellExperiment(Stemcell_meta, assay = "RNA")
  scoremodel_OLR <- psupertime(Species_OLR, Species_OLR$stage_order_sort, sel_genes = "list", gene_list = Commongene, n_folds = 5, method = "proportional")
  saveRDS(scoremodel_OLR, "Mouse_vRG_scoremodel.rds")
  # Perform linear regression for each gene
  expr_matrix <- t(scoremodel_OLR$x_data)
  stages <- scoremodel_OLR$proj_dt
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
  return(mean_expr)
}


Human_vRG_expr_mean <- get_mean_expr_after_OLR(Human_vRG, 60, set_stage_order = c("1","2","3","4","5"))
Human_ExN_IPC_expr_mean <- get_mean_expr_after_OLR(Human_ExN_IPC, 60, set_stage_order = c("1","2","3","4","5"))
Human_oRG_expr_mean <- get_mean_expr_after_OLR(Human_oRG, 60, set_stage_order = c("1","2","3","4","5"))

Macaque_vRG_expr_mean <- get_mean_expr_after_OLR(Macaque_vRG, 100, set_stage_order = c("1","2","3"))
Macaque_ExN_IPC_expr_mean <- get_mean_expr_after_OLR(Macaque_ExN_IPC, 80, set_stage_order = c("1","2","3"))
Macaque_oRG_expr_mean <- get_mean_expr_after_OLR(Macaque_oRG, 100, set_stage_order = c("1","2","3"))

Mouse_vRG_expr_mean <- get_mean_expr_after_OLR(Mouse_vRG, 60, set_stage_order = c("1","2","3","4","5"))
Mouse_ExN_IPC_expr_mean <- get_mean_expr_after_OLR(Mouse_ExN_IPC, 60, set_stage_order = c("1","2","3","4","5"))
# Only have 7 cells here
#Mouse_oRG_expr_mean <- get_mean_expr_after_OLR(Mouse_oRG, 10, set_stage_order = c("1","2","3","4","5"))

Ferret_vRG_expr_mean <- get_mean_expr_after_OLR(Ferret_vRG, 100, set_stage_order = c("1","2","3"))
Ferret_ExN_IPC_expr_mean <- get_mean_expr_after_OLR(Ferret_ExN_IPC, 100, set_stage_order = c("1","2","3"))
# Only have 161 cells
#Ferret_oRG_expr_mean <- get_mean_expr_after_OLR(Ferret_oRG, 100, set_stage_order = c("1","2","3"))

Rat_vRG_expr_mean <- get_mean_expr_after_OLR(Rat_vRG, 100, set_stage_order = c("1","2","3"))
Rat_ExN_IPC_expr_mean <- get_mean_expr_after_OLR(Rat_ExN_IPC, 100, set_stage_order = c("1","2","3"))
# Only have 365 cells
#Rat_oRG_expr_mean <- get_mean_expr_after_OLR(Rat_oRG, 100, set_stage_order = c("1","2","3"))


saveRDS(Human_vRG_expr_mean, "Human_vRG_expr_mean.rds")
saveRDS(Human_ExN_IPC_expr_mean, "Human_ExN_IPC_expr_mean.rds")
saveRDS(Human_oRG_expr_mean, "Human_oRG_expr_mean.rds")

saveRDS(Macaque_vRG_expr_mean, "Macaque_vRG_expr_mean.rds")
saveRDS(Macaque_ExN_IPC_expr_mean, "Macaque_ExN_IPC_expr_mean.rds")
saveRDS(Macaque_oRG_expr_mean, "Macaque_oRG_expr_mean.rds")

saveRDS(Ferret_vRG_expr_mean, "Ferret_vRG_expr_mean.rds")
saveRDS(Ferret_ExN_IPC_expr_mean, "Ferret_ExN_IPC_expr_mean.rds")

saveRDS(Mouse_vRG_expr_mean, "Mouse_vRG_expr_mean.rds")
saveRDS(Mouse_ExN_IPC_expr_mean, "Mouse_ExN_IPC_expr_mean.rds")

saveRDS(Rat_vRG_expr_mean, "Rat_vRG_expr_mean.rds")
saveRDS(Rat_ExN_IPC_expr_mean, "Rat_ExN_IPC_expr_mean.rds")
