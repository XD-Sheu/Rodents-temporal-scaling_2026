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
library(rgl)
library(plotly)


#Get the score model contain every orthologs
setwd("F:/Single cell analysis/Mouse and rat/Network and trend/objects")
Conserved_age_gene <- read_rds("shared_genes.rds")


#Read all orthologs containing score model_output from OLR model
Mouse_scoremodel <- read_rds("Mouse_vRG_scoremodel.rds")
Rat_scoremodel <- read_rds("Rat_vRG_scoremodel.rds")

plot_labels_over_psupertime_homemade(Mouse_scoremodel, palette = "YLGnBu")
plot_identified_genes_over_psupertime_homemade (Mouse_scoremodel, n_to_plot=20)

plot_labels_over_psupertime_homemade(Rat_scoremodel, palette = "YLGnBu")
plot_identified_genes_over_psupertime_homemade (Rat_scoremodel, n_to_plot=20)

scale_values <- function(x){(x-min(x))/(max(x)-min(x))}

expGlobalmouse <- t(Mouse_scoremodel$x_data)
rownames(expGlobalmouse) <- toupper(rownames(expGlobalmouse)) #uppcase mouse genes
trajmouse <- scale_values(as.numeric(t(Mouse_scoremodel$proj_dt$psuper)))
names(trajmouse) <- colnames(expGlobalmouse)

expGlobalrat <- t(Rat_scoremodel$x_data)
rownames(expGlobalrat) <- toupper(rownames(expGlobalrat)) #uppcase rat genes
trajrat <- scale_values(as.numeric(t(Rat_scoremodel$proj_dt$psuper)))
names(trajrat) <- colnames(expGlobalrat)


numPts = 500
winSize = 0.3

interGlobalmouse = cellAlign::interWeights(expDataBatch = expGlobalmouse, trajCond = trajmouse,
                                           winSz = winSize, numPts = numPts)

interGlobalrat = cellAlign::interWeights(expDataBatch = expGlobalrat, trajCond = trajrat,
                                           winSz = winSize, numPts = numPts)

saveRDS(interGlobalmouse, "interGlobalmouse_25_02_08.rds")
saveRDS(interGlobalrat, "interGlobalrat_25_02_08.rds")

sharedMarkers_mo_r = intersect(rownames(expGlobalmouse), rownames(expGlobalrat))


#scale the interpolated data (Recommended):
interScaledGlobalmouse = cellAlign::scaleInterpolate(interGlobalmouse)
interScaledGlobalrat = cellAlign::scaleInterpolate(interGlobalrat)


#####################################################################################
#perform global alignment

#Ref are the horizontal, the second passed when call the globalAlign function

#Mouse vs. Rat
rat_mouse = globalAlign(interScaledGlobalrat$scaledData[sharedMarkers_mo_r,], interScaledGlobalmouse$scaledData[sharedMarkers_mo_r,],
                        scores = list(query = interScaledGlobalrat$traj, ref = interScaledGlobalmouse$traj), sigCalc = T, numPerm = 500)
plotAlign(rat_mouse)
plotAlign_homemade(rat_mouse)
mouse_human_map = mapRealDataGlobal(rat_mouse,intTrajQuery = interScaledGlobalrat$traj, realTrajQuery = trajrat,
                                    intTrajRef = interScaledGlobalmouse$traj, realTrajRef = trajmouse)
plotMapping(mouse_human_map)


#####################################################################################
#perform local alignment of all shared genes:
Thresh=0.25
numPts = 200
alignment = localAlign(interScaledGlobalHuman$scaledData[sharedMarkers,],interScaledGlobalmacaque$scaledData[sharedMarkers,],threshPercent = Thresh, numPerm = 200)


#Mapping!
mapping = mapRealDataGlobal(alignment,intTrajQuery = interScaledGlobalmacaque$traj, realTrajQuery = trajmacaque,
                            intTrajRef = interScaledGlobalHuman$traj, realTrajRef = trajHuman)
plotMapping(mapping)


#Plot align
plotAlign_homemade <- function(alignment){
  costMat = alignment$localCostMatrix
  costMat = t(apply(costMat,1,function(x){return(as.numeric(x))}))
  linearInd = sub2ind(nrow(costMat), alignment$align[[1]]$index1, alignment$align[[1]]$index2)
  #costMat[linearInd] = NA
  costMat = data.frame(costMat, row.names=1:nrow(costMat))
  colnames(costMat) = 1:ncol(costMat)
  #for global alignment, where there is a pseudotime shift vector:
  if(!is.null(alignment$ptShift)){
    annotCols = data.frame(ptShift = abs(alignment$ptShift), sign = factor(sign(alignment$ptShift)),row.names = colnames(costMat))
    pheatmap(costMat, cluster_cols = F, cluster_rows=F, border_color = NA,
             main = 'alignment plot',
             show_rownames = F, show_colnames = F, annotation_col = annotCols, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
  }else{
    pheatmap(costMat, cluster_cols = F, cluster_rows=F, border_color = NA,
             main = 'alignment plot', show_rownames = F, show_colnames = F)
  }
  
  return(NA)
}


#align in a local way
