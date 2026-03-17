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
library(cellAlign)
library(dtw)
library(SingleCellExperiment)
library(scater)
library(tidyr)
library(destiny)
library(ggbeeswarm)
library(ggthemes)
library(reshape2)
library(pheatmap)
library(scales)
library(psupertime)



#SINGLE CELL DISTANCE MATRIX
data_integrated@active.assay <- "integrated"
RAWMATRIX <- as.data.frame(GetAssayData(subset(data_integrated, subset = Celltype %in% c("Neuroepithelium cell (NEC)","Radial glia cell (RGC)"))))
mousescore <- scoremodel$proj_dt #psuper project trained with all mouse stem cell include IP, n_folds = 10
ratscore <- scoremodel_v1$proj_dt #psuper project trained with all rat stem cell include IP, n_folds = 10
mouse_sorted <- mousescore[order(mousescore$psuper),]
rat_sorted <- ratscore[order(ratscore$psuper),]

distanceM <- dist(t(RAWMATRIX))
MM <- as.matrix(distanceM)  #can not directly transform dist object to dataframe. Needed to be firstly convert to matrix
heatmapM <- MM[mouse_sorted$cell_id, rat_sorted$cell_id]

col_fun = colorRamp2(seq(min(rat_sorted$psuper), max(rat_sorted$psuper), length.out = 30), hcl_palette = "Blues", reverse = T)
col_ha = HeatmapAnnotation(foo = rat_sorted$psuper,
                           bar = rat_sorted$label_input,
                           col = list(foo = col_fun,
                                      bar = c("14" = "#EFCDE3", "16" = "#D683B8", "18" = "#AE1877")))

col_fun = colorRamp2(seq(min(mouse_sorted$psuper), max(mouse_sorted$psuper), length.out = 30), hcl_palette = "Blues", reverse = T)
row_ha = rowAnnotation(foo = mouse_sorted$psuper,
                           bar = mouse_sorted$label_input,
                           col = list(foo = col_fun,
                                      bar = c("12" = "#EFCDE3", "14" = "#D683B8", "16" = "#AE1877")))

myBreaks <- seq(47, 81, length.out=100)#min(heatmapM) max(heatmapM)
hh <- Heatmap(as.matrix(heatmapM), cluster_rows = F, row_labels = rownames(heatmapM), show_row_names = F,row_names_side = c("left"), cluster_columns = F,
              column_labels = colnames(heatmapM), show_column_names = FALSE ,name = "Distance", clustering_distance_rows = "spearman", top_annotation = col_ha, left_annotation = row_ha,
              # a lot of decorations
              col=colorRamp2(myBreaks, hcl_palette = "Blue-Red 2", reverse = T), use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Gaussian")
hh

#IF group based on cluster is required (with their group names)
Heatmap(matrix(rnorm(100), 10), name = "mat",
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:4), labels = c("group1", "group2", "group3"), labels_gp = gpar(col = "white", fontsize = 10))), column_km = 3,
        left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 2:4), labels = c("group1", "group2", "group3"), labels_gp = gpar(col = "white", fontsize = 10))), row_km = 3)
#If want to mark severl gene in heatmap
m = matrix(rnorm(1000), nrow = 100)
rownames(m) = 1:100
ha = rowAnnotation(foo = anno_mark(at = c(1:4, 20, 60, 97:100), labels = month.name[1:10]))
Heatmap(m, name = "mat", cluster_rows = FALSE, right_annotation = ha, row_names_side = "left", row_names_gp = gpar(fontsize = 4))
Heatmap(m, name = "mat", cluster_rows = FALSE, right_annotation = ha, row_names_side = "left", row_names_gp = gpar(fontsize = 4), row_km = 4)
#GO TERM with HEATMAP
mat = matrix(rnorm(100*10), nrow = 100)
split = sample(letters[1:5], 100, replace = TRUE)
sentences = lapply(unique(split), function(x) {random_text(3, 8)})
names(sentences) = unique(split)
Heatmap(mat, name = "mat", row_split = split, right_annotation = rowAnnotation(textbox = anno_textbox(split, sentences,word_wrap = TRUE, add_new_line = TRUE)))



#PSUDOBULK
#ALL CELL
counts_DISTANCE_M_R <- as.matrix(Average_data_integrated$RNA)
metadata <- read.delim("Condition info.csv", header = TRUE, row.names = 1, sep = ',')
dds_DISTANCE_M_R <- DESeqDataSetFromMatrix(countData = round(counts_DISTANCE_M_R[which(rowSums(counts_DISTANCE_M_R)>0),]), colData = metadata[colnames(counts_DISTANCE_M_R),], design = ~Species)
ddsDE_DISTANCE_M_R <- DESeq(dds_DISTANCE_M_R)
vsdata_DISTANCE_M_R <- rlog(ddsDE_DISTANCE_M_R, blind=T)

ALLCELL_Dists <- as.matrix(dist(t(assay(vsdata_DISTANCE_M_R)), method = "euclidean"))
write.csv(ALLCELL_Dists, "ALLCELL distance.csv")
ALLCELL_pairwise <- read.delim("ALLCELL distance pairwise.csv", header = TRUE, row.names = 1, sep = ',')

mat.scaled <- t(apply(ALLCELL_pairwise, 1, scale)) 
colnames(mat.scaled) <- colnames(ALLCELL_pairwise)
colorcostoumed <- colorRampPalette(c("red", 'white' ,"blue"), alpha = F)(100)
myBreaks <- seq(-1, 1.35, length.out=100)
#raw matrix
Heatmap(mat.scaled, cluster_rows = F,
        row_labels = rownames(mat.scaled),
        row_names_side = c("left"),
        column_names_side = c("top"),
        row_names_gp = gpar(fontsize = 9),
        cluster_columns = F,
        column_labels = colnames(mat.scaled), name = "Distance",
        clustering_distance_rows = "pearson",
        # a lot of decorations
        col=colorRamp2(myBreaks, colorcostoumed),
        #col =  scale_color_viridis(discrete = TRUE, option = "D"),
        #column_dend_height=unit(20,"mm"),
        #column_names_gp = gpar(fontsize=12),
        # border_gp = gpar(col = "black", lty = 2),
        rect_gp = gpar(col = "white", lwd = 1),
        column_names_rot = 75,
)


#ONLY STEM CELL
Average_data_integrated_stem <- AverageExpression(subset(data_integrated, subset = Celltype %in% c("Neuroepithelium cell (NEC)","Cycling neuroepithelium cell (ccNEC)","Radial glia cell (RGC)","Cycling radial glia cell (ccRGC)",
                                                                                                   "Cycling intermediate progenitor (ccIP)", "Intermediate progenitor (IP)")), assays ="RNA", group.by = c("species", "day"), return.seurat = F)
counts_DISTANCE_M_R_stem <- as.matrix(Average_data_integrated_stem$RNA)
metadata <- read.delim("Condition info.csv", header = TRUE, row.names = 1, sep = ',')
dds_DISTANCE_M_R_stem <- DESeqDataSetFromMatrix(countData = round(counts_DISTANCE_M_R_stem[which(rowSums(counts_DISTANCE_M_R_stem)>0),]), colData = metadata[colnames(counts_DISTANCE_M_R_stem),], design = ~Species)
ddsDE_DISTANCE_M_R_stem <- DESeq(dds_DISTANCE_M_R_stem)
vsdata_DISTANCE_M_R_stem <- rlog(ddsDE_DISTANCE_M_R_stem, blind=T)

stem_Dists <- as.matrix(dist(t(assay(vsdata_DISTANCE_M_R_stem)), method = "euclidean"))
write.csv(stem_Dists, "STEM distance.csv")
stem_pairwise <- read.delim("STEM distance pairwise.csv", header = TRUE, row.names = 1, sep = ',')
mat.scaled <- t(apply(stem_pairwise, 1, scale)) 
colnames(mat.scaled) <- colnames(stem_pairwise)
colorcostoumed <- colorRampPalette(c("red", 'white' ,"blue"), alpha = F)(100)
myBreaks <- seq(-1, 1.35, length.out=100)
#raw matrix
Heatmap(mat.scaled, cluster_rows = F,
        row_labels = rownames(mat.scaled),
        row_names_side = c("left"),
        column_names_side = c("top"),
        row_names_gp = gpar(fontsize = 9),
        cluster_columns = F,
        column_labels = colnames(mat.scaled), name = "Distance",
        clustering_distance_rows = "pearson",
        # a lot of decorations
        col=colorRamp2(myBreaks, colorcostoumed),
        #col =  scale_color_viridis(discrete = TRUE, option = "D"),
        #column_dend_height=unit(20,"mm"),
        #column_names_gp = gpar(fontsize=12),
        # border_gp = gpar(col = "black", lty = 2),
        rect_gp = gpar(col = "white", lwd = 1),
        column_names_rot = 75,
)



