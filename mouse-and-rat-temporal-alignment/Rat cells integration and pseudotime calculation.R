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
library(ComplexHeatmap)
library(DESeq2)
library(clusterProfiler)
library(wordcloud)
library(enrichplot)
library(ggupset)
library(scplotter)
library(EnhancedVolcano)
library(RCurl)
library(AnnotationHub)
library(ensembldb)
library(SeuratWrappers)
library(patchwork)


Rat_cells <- subset(data_integrated, subset = species %in% "rat")
setwd("F:/Single cell analysis/Mouse and rat/Network and trend/objects")

#Score cell cycle phase
cc_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Mus_musculus.csv") 
cell_cycle_genes <- read.csv(text = cc_file)
ah <- AnnotationHub()
ahDb <- query(ah, pattern = c("Mus musculus", "EnsDb"), ignore.case = TRUE)
id <- ahDb %>% mcols() %>% rownames() %>% tail(n = 1)
edb <- ah[[id]]
annotations <- genes(edb, return.type = "data.frame")
annotations <- annotations %>% dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)
cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))
s_genes <- cell_cycle_markers %>% dplyr::filter(phase == "S") %>% pull("gene_name")
g2m_genes <- cell_cycle_markers %>% dplyr::filter(phase == "G2/M") %>% pull("gene_name")

Rat_cells <- CellCycleScoring(Rat_cells, s.features = s_genes, g2m.features = g2m_genes, set.ident = TRUE)
Rat_cells$CC.Difference <- Rat_cells$S.Score - Rat_cells$G2M.Score


Rat_cells@active.assay <- "RNA"
features <- SplitObject(Rat_cells, split.by = "day")
features <- lapply(X = features,
                   FUN = SCTransform,
                   method = "glmGamPoi",
                   vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt", "CC.Difference"),
                   return.only.var.genes = FALSE)
hvg_stemcell <- SelectIntegrationFeatures(features[grep("E14|E16|E18", names(features))], nfeatures = 4000) %>% unique()

for (i in seq_along(features)) {
  if ("integrated" %in% Assays(features[[i]])) {
    features[[i]][["integrated"]] <- NULL
  }
}

Rat_stemcell_regress <- merge(x = features[[1]], y = features[2:length(features)], merge.data=TRUE)
VariableFeatures(Rat_stemcell_regress) <- hvg_stemcell
Rat_stemcell_regress <- PrepSCTFindMarkers(Rat_stemcell_regress)
DefaultAssay(Rat_stemcell_regress) <- "SCT"

dayl <- c("1","2","3")
names(dayl) <- c("E14", "E16", "E18")
Idents(Rat_stemcell_regress) <- "day"
Rat_stemcell_regress <- RenameIdents(Rat_stemcell_regress, dayl)
Rat_stemcell_regress$stage_order <- Idents(Rat_stemcell_regress)

saveRDS(Rat_stemcell_regress, "Rat_stemcell_regress_25_08_06.rds")
Rat_stemcell_regress <- read_rds("Rat_stemcell_regress_25_08_06.rds")

hvg_Rat <- VariableFeatures(Rat_stemcell_regress)
write.csv(hvg_Rat, "hvg_Rat.csv")
write.csv(rownames(Rat_stemcell_regress), "genemodel_Rat.csv")

#integration and clustering
Rat_stemcell_regress <- RunPCA(Rat_stemcell_regress, verbose = FALSE, npcs = 50)
Rat_stemcell_regress <- FindNeighbors(Rat_stemcell_regress, dims = 1:50, reduction = "pca")
Rat_stemcell_regress <- FindClusters(Rat_stemcell_regress, resolution = 2.2, cluster.name = "unintegrated_clusters_stem")
Rat_stemcell_regress <- RunUMAP(Rat_stemcell_regress, reduction = "pca", dims = 1:30, reduction.name = "umap.unintegrated.stem")
stemcell_color <- c("Immune-and-blood" = "#f2f0e1",  "vRG" = "#ddc49e","oRG" = "#46a166","tRG" = "#92ec4a",  "Glia-IPC" = "#a6c4d1", "ExN-IPC" = "#e98074",  "OPC" = "#897dcb", "Astro-imma" = "#1e7894", "Astrocyte" = "#1e7894")

DimPlot(Rat_stemcell_regress, reduction = "umap.unintegrated.stem", group.by = "stemcell_type", cols = stemcell_color) & theme_void() 
stage_color <- hcl.colors(4, "Blues", rev = T)[2:4]
names(stage_color) <- c("E14","E16","E18")
DimPlot(Rat_stemcell_regress, reduction = "umap.unintegrated.stem", group.by = "day", cols = stage_color) & theme_void() 
FeaturePlot(Rat_stemcell_regress, features = Genename_stem, reduction = "umap.unintegrated.stem") & scale_colour_gradientn(colours = colorplate_green) & theme_void() & theme(legend.position="none") & theme(legend.position="none",text = element_text(size = 11))

#CCA 
Rat_stemcell_regress <- IntegrateLayers(object = NormalizeData(Rat_stemcell_regress), method = CCAIntegration, normalization.method = "SCT",orig.reduction = "pca", new.reduction = "integrated.cca.stem", verbose = T)
Rat_stemcell_regress <- FindNeighbors(Rat_stemcell_regress, reduction = "integrated.cca.stem", dims = 1:50)
Rat_stemcell_regress <- FindClusters(Rat_stemcell_regress, resolution = 1.2, cluster.name = "cca_clusters_stem")
Rat_stemcell_regress <- RunUMAP(Rat_stemcell_regress, reduction = "integrated.cca.stem", dims = 1:20, reduction.name = "umap.cca.stem")
#RPCA 
Rat_stemcell_regress <- IntegrateLayers(object = Rat_stemcell_regress, method = RPCAIntegration, normalization.method = "SCT", orig.reduction = "pca", new.reduction = "integrated.rpca.stem", verbose = T)
Rat_stemcell_regress <- FindNeighbors(Rat_stemcell_regress, reduction = "integrated.rpca.stem", dims = 1:50)
Rat_stemcell_regress <- FindClusters(Rat_stemcell_regress, resolution = 1.2, cluster.name = "rpca_clusters_stem")
Rat_stemcell_regress <- RunUMAP(Rat_stemcell_regress, reduction = "integrated.rpca.stem", dims = 1:20, reduction.name = "umap.rpca.stem")
#HARMONY
Rat_stemcell_regress <- IntegrateLayers(object = Rat_stemcell_regress, method = HarmonyIntegration, assay = "SCT", group.by.vars = "day", orig.reduction = "pca", new.reduction = "harmony.stem", verbose = T)
Rat_stemcell_regress <- FindNeighbors(Rat_stemcell_regress, reduction = "harmony.stem", dims = 1:50)
Rat_stemcell_regress <- FindClusters(Rat_stemcell_regress, resolution = 2, cluster.name = "harmony_clusters_stem")
Rat_stemcell_regress <- RunUMAP(Rat_stemcell_regress, reduction = "harmony.stem", dims = 1:30, reduction.name = "umap.harmony.stem")
DimPlot(Rat_stemcell_regress, reduction = "umap.harmony.stem", group.by = "day", cols = stage_color) & theme_void() 
FeaturePlot(Rat_stemcell_regress, features = Genename_stem, reduction = "umap.harmony.stem") & scale_colour_gradientn(colours = colorplate_green) & theme_void() & theme(legend.position="none") & theme(legend.position="none",text = element_text(size = 11))
#MNN
DefaultAssay(Rat_stemcell_regress) <- "SCT"
Rat_stemcell_regress <- IntegrateLayers(object = Rat_stemcell_regress, method = FastMNNIntegration, new.reduction = "integrated.mnn.stem", batch = Rat_stemcell_regress$day, verbose = T)
Rat_stemcell_regress <- FindNeighbors(Rat_stemcell_regress, reduction = "integrated.mnn.stem", dims = 1:50)
Rat_stemcell_regress <- FindClusters(Rat_stemcell_regress, resolution = 2, cluster.name = "mnn_clusters_stem")
Rat_stemcell_regress <- RunUMAP(Rat_stemcell_regress, reduction = "integrated.mnn.stem", dims = 1:20, reduction.name = "umap.mnn.stem")

saveRDS(Rat_stemcell_regress, "Rat_stemcell_regress_25_02_26.rds")

stage_color <- hcl.colors(11, "Blues", rev = T)
names(stage_color) <- c("E10","E11","E12","E13","E14","E15","E16","E17","E18","P1","P4")
DefaultAssay(Rat_stemcell_regress) <- "SCT"
P1 <- DimPlot(Rat_stemcell_regress, reduction = "umap.unintegrated.stem", group.by = "day", cols = stage_color, pt.size = 0.8) & theme_void() & theme(legend.position="none")
P2 <- DimPlot(Rat_stemcell_regress, reduction = "umap.cca.stem", group.by = "day", cols = stage_color, pt.size = 0.8) & theme_void() & theme(legend.position="none")
P3 <- DimPlot(Rat_stemcell_regress, reduction = "umap.rpca.stem", group.by = "day", cols = stage_color, pt.size = 0.8) & theme_void() & theme(legend.position="none")
P4 <- DimPlot(Rat_stemcell_regress, reduction = "umap.harmony.stem", group.by = "day", cols = stage_color, pt.size = 0.8) & theme_void() & theme(legend.position="none")
P5 <- DimPlot(Rat_stemcell_regress, reduction = "umap.mnn.stem", group.by = "day", cols = stage_color, pt.size = 0.8) & theme_void() & theme(legend.position="none")

P1 | P2 | P3 | P4 | P5

colorplate_green <- c("#cdd5d5","#FCFFDD", "#77D1B5", "#00A6AE", "#006AA8", "#26185F")
colorplate_red <- c("#cdd5d5","#FCF5F2", "#FFCDC8", "#F893AE", "#DF44A8", "#95148D", "#490062")

Genename_stem <- c("Hmga2","Pax6","Sox2","Mki67","Ccnb2","Ccnd1","Hopx","Eomes","Aldh1l1","Apoe","Olig2","Pdgfra")
Rat_stemcell_regress@active.assay <- "SCT"
FeaturePlot(Rat_stemcell_regress, features = Genename_stem, reduction = "umap.unintegrated.stem") & scale_colour_gradientn(colours = colorplate_green) & theme_void() & theme(legend.position="none")

DimPlot(Rat_stemcell_regress, reduction = "umap.unintegrated.stem", group.by = c("S.Score","G2M.Score"))  & theme_void()



#Feature genes in different reductions
colorplate_green <- c("#cdd5d5","#FCFFDD", "#77D1B5", "#00A6AE", "#006AA8", "#26185F")
colorplate_red <- c("#cdd5d5","#FCF5F2", "#FFCDC8", "#F893AE", "#DF44A8", "#95148D", "#490062")

Genename_stem <- c("Hmga2","Hes1","Eomes","Dcx","Hopx","Tnc","Lifr","Cryab","Satb2","Tbr1","Bcl11b","Cux1","Egfr","Gfap","Olig2","Pdgfra")

Rat_stemcell_regress@active.assay <- "SCT"
FeaturePlot(Rat_stemcell_regress, features = Genename_stem, reduction = "umap.unintegrated.stem") & scale_colour_gradientn(colours = colorplate_green) & theme_void() & theme(legend.position="none") & theme(legend.position="none",text = element_text(size = 11))
FeaturePlot(Rat_stemcell_regress, features = Genename_stem, reduction = "umap.cca.stem") & scale_colour_gradientn(colours = colorplate_red) & theme_void() & theme(legend.position="none") & theme(legend.position="none",text = element_text(size = 11))
FeaturePlot(Rat_stemcell_regress, features = Genename_stem, reduction = "umap.rpca.stem") & scale_colour_gradientn(colours = colorplate_green) & theme_void() & theme(legend.position="none") & theme(legend.position="none",text = element_text(size = 11))
FeaturePlot(Rat_stemcell_regress, features = Genename_stem, reduction = "umap.harmony.stem") & scale_colour_gradientn(colours = colorplate_red) & theme_void() & theme(legend.position="none") & theme(legend.position="none",text = element_text(size = 11))
FeaturePlot(Rat_stemcell_regress, features = Genename_stem, reduction = "umap.mnn.stem") & scale_colour_gradientn(colours = colorplate_green) & theme_void() & theme(legend.position="none",text = element_text(size = 11))



#Determine clusters
Idents(Rat_stemcell_regress) <- "unintegrated_clusters_stem"
DEG_markers <- FindAllMarkers(Rat_stemcell_regress, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.2,test.use = "wilcox")
DEG_markers %>% group_by(cluster) %>% dplyr::filter(p_val_adj < 0.05) %>% dplyr::filter(pct.1 > 0.3) %>% top_n(n = 50, wt = avg_log2FC)  %>% ungroup() -> top100_DEG_markers_Rat_celltype
DotPlot(Rat_stemcell_regress, features = Genename_stem, cols = hcl.colors(2, palette = "Berlin"), dot.scale = 5, scale.by = "radius", cluster.idents = T) + RotatedAxis() + theme(axis.title.y=element_blank(), axis.title.x=element_blank()) 
DimPlot(Rat_stemcell_regress, reduction = "umap.unintegrated.stem", group.by = "unintegrated_clusters_stem", cols = hcl.colors(35, palette = "Dark3", rev = T), label = T, label.color = "white", label.box = T, repel = T) & theme_void()
setwd("F:/Single cell analysis/Mouse and rat/Network and trend/objects")
write.csv(top100_DEG_markers_Rat_celltype, "top100_DEG_markers_Rat.csv")


#unintegrated cell cluster
Idents(Rat_stemcell_regress) <- "unintegrated_clusters_stem"
setwd("F:/Single cell analysis/Mouse and rat/Network and trend/objects")
new.cluster.ids <- read.csv("Cell type annotation_Rat_unintegrated.csv", header = F)$V2
names(new.cluster.ids) <- as.character(read.csv("Cell type annotation_Rat_unintegrated.csv", header = F)$V1)
Rat_stemcell_regress <- RenameIdents(Rat_stemcell_regress, new.cluster.ids)
Rat_stemcell_regress$stemcell_type <- Idents(Rat_stemcell_regress)
Idents(Rat_stemcell_regress) <- "stemcell_type"
#Remove the DCX soly expressed differentiated newborn neuron
Rat_stemcell_regress <- subset(Rat_stemcell_regress, 
                                 subset = unintegrated_clusters_stem %in% c(50, 51, 64, 73, 89, 77, 54, 87, 57), 
                                 invert = TRUE)

DimPlot(Rat_stemcell_regress, reduction = "umap.mnn.stem", group.by = "stemcell_type", cols = stemcell_color) & theme_void() 


#Validation of clustering
Rat_stemcell_regress@active.assay <- "SCT"
Genename_stem <- c("Ran","Hmga2","Dach1","Lef1","Eomes","Mllt3","Auts2","Myt1l","Tnc","Notch2","Hopx","Ccdc175","Cryab","Fstl1","Adgrv1","Prr16","Egfr","Spry2","Map3k1","Dpf3","Aqp4","Gfap","Npnt","Rbms3","Pdgfra","Pcdh15","Car10","Dock10")
stemcell_color <- c("vRG" = "#ddc49e","ExN-IPC" = "#e98074",  "oRG" = "#46a166","tRG" = "#82d83d",  "Glia-IPC" = "#a6c4d1", "Astro-imma" = "#1e7894", "OPC" = "#897dcb")

Idents(Rat_stemcell_regress) <- "stemcell_type"
stemcell_type_sorted <- c("vRG","ExN-IPC",  "oRG", "tRG", "Glia-IPC", "Astro-imma", "OPC")
Idents(Rat_stemcell_regress) <- factor(Idents(Rat_stemcell_regress), levels= stemcell_type_sorted)
Rat_stemcell_regress$stemcell_type_sorted <- Idents(Rat_stemcell_regress)
Idents(Rat_stemcell_regress) <- "stemcell_type_sorted"
Rat_stemcell_regress$Age_numeric <- sapply(Rat_stemcell_regress$day, function(x) {
  num <- as.numeric(gsub("[EP]", "", x))
  if (grepl("^P", x)) {
    return(num + 19)
  } else {
    return(num)
  }
})

ht <- GroupHeatmap(
  assay = "SCT",
  slot  =  "data",
  srt = Rat_stemcell_regress,
  flip = T,
  features = Genename_stem,
  group.by = c("stemcell_type_sorted"),
  group_palcolor = list(stemcell_color),
  fill_palcolor = list(c("#ddc49e","#ddc49e","#ddc49e","#ddc49e","#e98074","#e98074","#e98074","#e98074","#46a166","#46a166","#46a166","#46a166",
                         "#82d83d","#82d83d","#82d83d","#82d83d","#a6c4d1","#a6c4d1","#a6c4d1","#a6c4d1","#1e7894","#1e7894","#1e7894","#1e7894","#897dcb","#897dcb","#897dcb","#897dcb")),
  cell_annotation = c("Age_numeric", "Pcna", "Phase"),
  cell_annotation_palcolor = list(stemcell_color, stemcell_color, c("#9a9a80","#e6e6df","#353500")),
  row_names_side = "left",
  add_violin = TRUE, 
  #add_dot = TRUE, add_reticle = TRUE,
  exp_cutoff = 0,
  cluster_columns = F, show_column_names = F, show_row_names = T
)
print(ht$plot)



ht <- GroupHeatmap(
  assay = "SCT",
  slot = "data",
  srt = Rat_stemcell_regress,
  cluster_rows = T,
  features = Genename_stem,
  group.by = c("unintegrated_clusters_stem"),
  cell_annotation = c("Age_numeric"),
  row_names_side = "left",
  show_column_names = T
)
print(ht$plot)






#Psupertime based trajecotry analysis
#change stage info to ordinal sequence
stage_order <-  c("1","2","3","4","5","6","7","8","9","10","11")
names(stage_order) <- c("E10","E11","E12","E13","E14","E15","E16","E17","E18","P1","P4")
Idents(Rat_stemcell_regress) <- "stage"
Rat_stemcell_regress <- RenameIdents(Rat_stemcell_regress, stage_order)
Rat_stemcell_regress$stage_order <- Idents(Rat_stemcell_regress)

Rat_stemcell_regress$metacell <- NA

stages <- unique(Rat_stemcell_regress$stage_order)

for (stage_i in stages) {
  message("Processing stage: ", stage_i)
  seurat_subset <- subset(Rat_stemcell_regress, subset = stage_order == stage_i)
  seurat_subset <- FindClusters(seurat_subset, resolution = 2, verbose = FALSE)
  
  initial_clusters <- seurat_subset$seurat_clusters
  metacells <- split(colnames(seurat_subset), initial_clusters)
  message("Initial FindClusters produced ", length(metacells), " clusters")
  
  while (length(metacells) < 50) {
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
  metacell_labels <- paste0(stage_i, "_", 1:length(metacells))
  for (i in 1:length(metacells)) {
    Rat_stemcell_regress$metacell[metacells[[i]]] <- metacell_labels[i]
  }
  message("Stage ", stage_i, " completed with ", length(metacells), " metacells")
}

averaged_exp <- AggregateExpression(Rat_stemcell_regress, group.by = "metacell")
averaged_exp <- averaged_exp[[1]]
colnames(averaged_exp) <- sub("^g(.*)", "\\1", colnames(averaged_exp))
colnames(averaged_exp) <- gsub("-", "_", colnames(averaged_exp))

metacell_types <- sapply(unique(Rat_stemcell_regress$metacell), function(m) {
  cells <- which(Rat_stemcell_regress$metacell == m)
  cell_types <- Rat_stemcell_regress$stemcell_type[cells]
  tbl <- table(cell_types)
  names(tbl)[which.max(tbl)]
})

stages <- sapply(strsplit(colnames(averaged_exp), "_"), function(x) x[1])
meta_data <- data.frame(
  stage_order = stages,
  stemcell_type = metacell_types[colnames(averaged_exp)],
  row.names = colnames(averaged_exp)
)

Rat_stemcell_meta <- CreateSeuratObject(counts = averaged_exp, meta.data = meta_data, project = "Rat_stemcell_meta")
Rat_stemcell_meta <- NormalizeData(Rat_stemcell_meta, normalization.method = "LogNormalize", scale.factor = 10000)
DefaultAssay(Rat_stemcell_meta) <- "RNA"
Idents(Rat_stemcell_meta) <- "stage_order"
Rat_stemcell_meta <- FindVariableFeatures(Rat_stemcell_meta, selection.method = "vst", nfeatures = 2000)
hvg_stemcell <- unique(c(VariableFeatures(Rat_stemcell_regress),VariableFeatures(Rat_stemcell_regress)))
setwd("F:/Single cell analysis/Parallel development/Rat/Objects")
saveRDS(Rat_stemcell_meta,"Rat_stemcell_meta_25_05_20.rds")



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
setwd("F:/Single cell analysis/Parallel development/Paper figures/Figure 7/objects")

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
setwd("F:/Single cell analysis/Parallel development/Paper figures/Figure 2/OLR_pseudotime")
write.csv(top_500_genes, "Top 500 linear genes Rat.csv")
write.csv(results, "Linear genes all Rat.csv")
top_10_increasing <- top_500_genes[top_500_genes$coefficient > 0, ][1:40, ]$gene
top_10_decreasing <- top_500_genes[top_500_genes$coefficient < 0, ][1:40, ]$gene
Interest_20 <- c(top_10_decreasing, top_10_increasing)
Interest_20


#OLR model
Idents(Rat_stemcell_meta) <- "stage_order"
stage_order_sorted <- c("1","2","3","4","5","6","7","8","9","10","11")
Idents(Rat_stemcell_meta) <- factor(Idents(Rat_stemcell_meta), levels= stage_order_sorted)
Rat_stemcell_meta$stage_order_sorted <- Idents(Rat_stemcell_meta)
Idents(Rat_stemcell_meta) <- "stage_order_sorted"

Rat_OLR <- as.SingleCellExperiment(Rat_stemcell_meta, assay = "RNA")
setwd("F:/Single cell analysis/Parallel development/Rat/Objects")
conserved_gene_list <- read.csv("Atlas conserved age genes.csv")
conserved_gene_list <- conserved_gene_list$Rat.gene.symbol

scoremodel_Rat <- psupertime(Rat_OLR, Rat_OLR$stage_order_sorted, sel_genes = "list", gene_list = conserved_gene_list, n_folds = 3, method = "proportional")
#scoremodel_Rat <- psupertime(Rat_OLR, Rat_OLR$stage_order_sort, n_folds = 5, method = "proportional")

setwd("F:/Single cell analysis/Parallel development/Rat/Objects")
saveRDS(scoremodel_Rat,"scoremodel_Rat_25_05_20.rds")

plot_labels_over_psupertime_homemade(scoremodel_Rat, palette = "YLGnBu")
plot_identified_genes_over_psupertime_homemade (scoremodel_Rat, n_to_plot=63)
plot_specified_genes_over_psupertime_homemade(scoremodel_Rat, Interest_20, palette = "Blues")
plot_specified_genes_over_psupertime_homemade(scoremodel_Rat, c("Hmga2"), palette = "Blues")
view(scoremodel_Rat$beta_dt)
view(goanalysis$clusters_dt)
view(scoremodel$x_data)

setwd("F:/Single cell analysis/Parallel development/Rat/Objects")
scoremodel_Rat <- read_rds("scoremodel_Rat_25_05_20.rds")













