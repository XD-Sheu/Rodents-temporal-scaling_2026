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
Conserved_age_gene <- read_rds("Five_animal_shared_genes.rds")


#Read all orthologs containing score model
Human_scoremodel <- read_rds("Human_vRG_scoremodel.rds")
Macaque_scoremodel <- read_rds("Macaque_vRG_scoremodel.rds")
Ferret_scoremodel <- read_rds("Ferret_vRG_scoremodel.rds")
Mouse_scoremodel <- read_rds("Mouse_vRG_scoremodel.rds")
#Mouse_scoremodel <- read_rds("Mouse_vRG_scoremodel_alternative.rds")
Rat_scoremodel <- read_rds("Rat_vRG_scoremodel.rds")

scale_values <- function(x){(x-min(x))/(max(x)-min(x))}

expGlobalHuman <- t(Human_scoremodel$x_data)
trajHuman <- scale_values(as.numeric(t(Human_scoremodel$proj_dt$psuper)))
names(trajHuman) <- colnames(expGlobalHuman)

expGlobalmacaque <- t(Macaque_scoremodel$x_data)
trajmacaque <- scale_values(as.numeric(t(Macaque_scoremodel$proj_dt$psuper)))
names(trajmacaque) <- colnames(expGlobalmacaque)

expGlobalferret <- t(Ferret_scoremodel$x_data)
rownames(expGlobalferret) <- sub("\\..*", "", rownames(expGlobalferret)) #remove .BLAST in ferret genes
trajferret <- scale_values(as.numeric(t(Ferret_scoremodel$proj_dt$psuper)))
names(trajferret) <- colnames(expGlobalferret)

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
interGlobalHuman = cellAlign::interWeights(expDataBatch = expGlobalHuman, trajCond = trajHuman,
                                           winSz = winSize, numPts = numPts)

interGlobalmacaque = cellAlign::interWeights(expDataBatch = expGlobalmacaque, trajCond = trajmacaque,
                                             winSz = winSize, numPts = numPts)

interGlobalferret = cellAlign::interWeights(expDataBatch = expGlobalferret, trajCond = trajferret,
                                            winSz = winSize, numPts = numPts)

interGlobalmouse = cellAlign::interWeights(expDataBatch = expGlobalmouse, trajCond = trajmouse,
                                           winSz = winSize, numPts = numPts)

interGlobalrat = cellAlign::interWeights(expDataBatch = expGlobalrat, trajCond = trajrat,
                                           winSz = winSize, numPts = numPts)


#scale the interpolated data (Recommended):
interScaledGlobalHuman = cellAlign::scaleInterpolate(interGlobalHuman)
interScaledGlobalmacaque = cellAlign::scaleInterpolate(interGlobalmacaque)
interScaledGlobalferret = cellAlign::scaleInterpolate(interGlobalferret)
interScaledGlobalmouse = cellAlign::scaleInterpolate(interGlobalmouse)
interScaledGlobalrat = cellAlign::scaleInterpolate(interGlobalrat)


#####################################################################################
#perform global alignment

#Ref are the horizontal, the second passed when call the globalAlign function

#Human vs. macaque
macaque_human = globalAlign(interScaledGlobalmacaque$scaledData[sharedMarkers_all,], interScaledGlobalHuman$scaledData[sharedMarkers_all,],
                            scores = list(query = interScaledGlobalmacaque$traj, ref = interScaledGlobalHuman$traj), sigCalc = T, numPerm = 500)
plotAlign(macaque_human)
plotAlign_homemade(macaque_human)
macaque_human_map = mapRealDataGlobal(macaque_human,intTrajQuery = interScaledGlobalmacaque$traj, realTrajQuery = trajmacaque,
                                      intTrajRef = interScaledGlobalHuman$traj, realTrajRef = trajHuman)
plotMapping(macaque_human_map)


#Human vs. ferret
ferret_human = globalAlign(interScaledGlobalferret$scaledData[sharedMarkers_all,], interScaledGlobalHuman$scaledData[sharedMarkers_all,],
                           scores = list(query = interScaledGlobalferret$traj, ref = interScaledGlobalHuman$traj), sigCalc = T, numPerm = 500)
plotAlign(ferret_human)
plotAlign_homemade(ferret_human)
ferret_human_map = mapRealDataGlobal(ferret_human,intTrajQuery = interScaledGlobalferret$traj, realTrajQuery = trajferret,
                                     intTrajRef = interScaledGlobalHuman$traj, realTrajRef = trajHuman)
plotMapping(ferret_human_map)


#Human vs. mouse
mouse_human = globalAlign(interScaledGlobalmouse$scaledData[sharedMarkers_all,], interScaledGlobalHuman$scaledData[sharedMarkers_all,],
                          scores = list(query = interScaledGlobalmouse$traj, ref = interScaledGlobalHuman$traj), sigCalc = T, numPerm = 500)
plotAlign(mouse_human)
plotAlign_homemade(mouse_human)
mouse_human_map = mapRealDataGlobal(mouse_human,intTrajQuery = interScaledGlobalmouse$traj, realTrajQuery = trajmouse,
                                    intTrajRef = interScaledGlobalHuman$traj, realTrajRef = trajHuman)
plotMapping(mouse_human_map)

#Human vs. rat
rat_human = globalAlign(interScaledGlobalrat$scaledData[sharedMarkers_all,], interScaledGlobalHuman$scaledData[sharedMarkers_all,],
                          scores = list(query = interScaledGlobalrat$traj, ref = interScaledGlobalHuman$traj), sigCalc = T, numPerm = 500)
plotAlign(rat_human)
plotAlign_homemade(rat_human)
mouse_human_map = mapRealDataGlobal(rat_human,intTrajQuery = interScaledGlobalrat$traj, realTrajQuery = trajrat,
                                    intTrajRef = interScaledGlobalHuman$traj, realTrajRef = trajHuman)
plotMapping(mouse_human_map)

#Macaque vs. Ferret
ferret_macaque = globalAlign(interScaledGlobalferret$scaledData[sharedMarkers_all,], interScaledGlobalmacaque$scaledData[sharedMarkers_all,],
                             scores = list(query = interScaledGlobalferret$traj, ref = interScaledGlobalmacaque$traj), sigCalc = T, numPerm = 500)
plotAlign(ferret_macaque)
plotAlign_homemade(ferret_macaque)
ferret_macaque_map = mapRealDataGlobal(ferret_macaque,intTrajQuery = interScaledGlobalferret$traj, realTrajQuery = trajferret,
                                       intTrajRef = interScaledGlobalmacaque$traj, realTrajRef = trajmacaque)
plotMapping(ferret_macaque_map)


#Macaque vs. mouse
mouse_macaque = globalAlign(interScaledGlobalmouse$scaledData[sharedMarkers_all,], interScaledGlobalmacaque$scaledData[sharedMarkers_all,],
                            scores = list(query = interScaledGlobalmouse$traj, ref = interScaledGlobalmacaque$traj), sigCalc = T, numPerm = 500)
plotAlign(mouse_macaque)
plotAlign_homemade(mouse_macaque)
mouse_macaque_map = mapRealDataGlobal(mouse_macaque,intTrajQuery = interScaledGlobalmouse$traj, realTrajQuery = trajmouse,
                                      intTrajRef = interScaledGlobalmacaque$traj, realTrajRef = trajmacaque)
plotMapping(mouse_macaque_map)


#Macaque vs. rat
rat_macaque = globalAlign(interScaledGlobalrat$scaledData[sharedMarkers_all,], interScaledGlobalmacaque$scaledData[sharedMarkers_all,],
                        scores = list(query = interScaledGlobalrat$traj, ref = interScaledGlobalmacaque$traj), sigCalc = T, numPerm = 500)
plotAlign(rat_macaque)
plotAlign_homemade(rat_macaque)
mouse_human_map = mapRealDataGlobal(rat_macaque,intTrajQuery = interScaledGlobalrat$traj, realTrajQuery = trajrat,
                                    intTrajRef = interScaledGlobalmacaque$traj, realTrajRef = trajHuman)
plotMapping(mouse_human_map)


#Ferret vs. mouse
mouse_ferret = globalAlign(interScaledGlobalmouse$scaledData[sharedMarkers_all,], interScaledGlobalferret$scaledData[sharedMarkers_all,],
                           scores = list(query = interScaledGlobalmouse$traj, ref = interScaledGlobalferret$traj), sigCalc = T, numPerm = 500)
plotAlign(mouse_ferret)
plotAlign_homemade(mouse_ferret)
mouse_ferret_map = mapRealDataGlobal(mouse_ferret,intTrajQuery = interScaledGlobalmouse$traj, realTrajQuery = trajmouse,
                                     intTrajRef = interScaledGlobalferret$traj, realTrajRef = trajferret)
plotMapping(mouse_ferret_map)


#Ferret vs. rat
rat_ferret = globalAlign(interScaledGlobalrat$scaledData[sharedMarkers_all,], interScaledGlobalferret$scaledData[sharedMarkers_all,],
                        scores = list(query = interScaledGlobalrat$traj, ref = interScaledGlobalferret$traj), sigCalc = T, numPerm = 500)
plotAlign(rat_ferret)
plotAlign_homemade(rat_ferret)
mouse_human_map = mapRealDataGlobal(rat_ferret,intTrajQuery = interScaledGlobalrat$traj, realTrajQuery = trajrat,
                                    intTrajRef = interScaledGlobalferret$traj, realTrajRef = trajHuman)
plotMapping(mouse_human_map)


#Mouse vs. Rat
rat_mouse = globalAlign(interScaledGlobalrat$scaledData[sharedMarkers_all,], interScaledGlobalmouse$scaledData[sharedMarkers_all,],
                         scores = list(query = interScaledGlobalrat$traj, ref = interScaledGlobalmouse$traj), sigCalc = T, numPerm = 500)
plotAlign(rat_mouse)
plotAlign_homemade(rat_mouse)
mouse_human_map = mapRealDataGlobal(rat_mouse,intTrajQuery = interScaledGlobalrat$traj, realTrajQuery = trajrat,
                                    intTrajRef = interScaledGlobalmouse$traj, realTrajRef = trajmouse)
plotMapping(mouse_human_map)



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
