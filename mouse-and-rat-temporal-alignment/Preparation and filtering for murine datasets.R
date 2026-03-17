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

#Prepare for rat cells
rat_file <- read.delim(file = "F:/Single cell analysis/Mouse and rat/ExpMatrix_rat_MS446_mrg_RECODE_log2ss100k.csv",  
                       sep = ',', header = TRUE, row.names = 1, stringsAsFactors = FALSE)
rat <- CreateSeuratObject(counts = t(rat_file), project = "sc_rat_all_stages", min.cells = 0, min.features = 200)
rat$species <- "rat"
rat$day <- as.character(rat$orig.ident) #make the day information from factor to charcter type
rat[["percent.mt"]] <- PercentageFeatureSet(rat, pattern = "^mt-")
rat[["percent.ribo"]] <- PercentageFeatureSet(rat, pattern = "^Rp[sl]")

rat$log10GenesPerUMI <- log10(rat$nFeature_RNA) / log10(rat$nCount_RNA)
metadata_rat <- rat@meta.data 
metadata_rat %>% 
  ggplot(aes(x=day, fill=day)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

metadata_rat %>% 
  ggplot(aes(color=day, x=nCount_RNA, fill= day)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 15000) #cutoff

metadata_rat %>% 
  ggplot(aes(color=day, x=nFeature_RNA, fill= day)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 15000)

metadata_rat %>% 
  ggplot(aes(x=day, y=log10(nFeature_RNA), fill=day)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs nFeature_RNAs")

metadata_rat %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 15000) +
  geom_hline(yintercept = 15000) +
  facet_wrap(~day, nrow = 3)

metadata_rat %>% 
  ggplot(aes(color=day, x=percent.mt, fill=day)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 25)

metadata_rat %>%
  ggplot(aes(x=log10GenesPerUMI, color = day, fill=day)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.9)


VlnPlot(rat, features = c("nFeature_RNA","nCount_RNA"), pt.size = 0.1, cols = c("#67A0D9","#A0D967","#D967A0"), alpha = 0.05) & 
  theme(plot.title = element_text(size=10))


#filter rat cells
filtered_rat <- subset(x = rat, 
                        subset= (nCount_RNA >= 15000) & 
                          (nFeature_RNA >= 15000) & 
                          (log10GenesPerUMI > 0.90) & 
                          (percent.mt < 25) )










#mouse data
files <- list.files("F:/Single cell analysis/Mouse and rat/Mouse H5 data stage picked up")
input_files <- function(h5_path){
  dataname <- sub(".h5", "\\1", h5_path)
  day <- substr(dataname, 2, 3)
  data <- Read10X_h5(paste0("F:/Single cell analysis/Mouse and rat/Mouse H5 data stage picked up/", h5_path))
  
  data <- CreateSeuratObject(data, min.cells = 0, min.features = 200)
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")
  data[["percent.ribo"]] <- PercentageFeatureSet(data, pattern = "^Rp[sl]")
  data$day <- paste0("E", day)
  data$species <- "mouse"
  
  return(data)
}

mouse_list <- sapply(files, input_files)
#simple way to just merge
mouse <- merge(mouse_list[1]$E10.h5, y = mouse_list[2:length(mouse_list)])

mouse$log10GenesPerUMI <- log10(mouse$nFeature_RNA) / log10(mouse$nCount_RNA)
metadata_mouse <- mouse@meta.data 
metadata_mouse %>% 
  ggplot(aes(x=day, fill=day)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

metadata_mouse %>% 
  ggplot(aes(color=day, x=nCount_RNA, fill= day)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 1000) #cutoff

metadata_mouse %>% 
  ggplot(aes(color=day, x=nFeature_RNA, fill= day)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 800)

metadata_mouse %>% 
  ggplot(aes(x=day, y=log10(nFeature_RNA), fill=day)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs nFeature_RNAs")

metadata_mouse %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 1000) +
  geom_hline(yintercept = 1000) +
  facet_wrap(~day)

metadata_mouse %>% 
  ggplot(aes(color=day, x=percent.mt, fill=day)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 20)

metadata_mouse %>%
  ggplot(aes(x=log10GenesPerUMI, color = day, fill=day)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)


VlnPlot(mouse, features = c("nFeature_RNA","nCount_RNA"), pt.size = 0.1, cols = c("#67A0D9","#A0D967","#D967A0"), alpha = 0.05) & 
  theme(plot.title = element_text(size=10))

#Filtering mouse cells
filtered_mouse <- subset(x = mouse, 
                         subset= (nCount_RNA >= 1000) & 
                         (nFeature_RNA >= 800) & 
                         (log10GenesPerUMI > 0.75) & 
                         (percent.mt < 20))
filtered_mouse <- JoinLayers(filtered_mouse)







#merge two project
#simple way
#It seems they have different total read or read depths differs. Therefore use we use different scale.factor
filtered_mouse <- NormalizeData(object = filtered_mouse, normalization.method = "LogNormalize", scale.factor = 58900)
filtered_rat <- NormalizeData(object = filtered_rat, normalization.method = "LogNormalize", scale.factor = 58900)
common.features <- intersect(rownames(filtered_mouse), rownames(filtered_mouse))
data <- merge(filtered_mouse[common.features, ], filtered_rat[common.features, ], merge.data = T)
data <- JoinLayers(data)
VlnPlot(data, features = c("Gapdh","Ubc", "Actb", "Tubb5","Ywhaz"), group.by = "species", split.by = "day", pt.size = 0)



#Hardway
# Select the most variable features to use for integration
options(future.globals.maxSize = 8000 * 1024^2)
split_seurat <- SplitObject(merge(filtered_mouse, filtered_rat, merge.data = T), split.by = "species")
split_seurat <- split_seurat[c("mouse","rat")]

for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
  split_seurat[[i]] <- SCTransform(split_seurat[[i]])
}

integ_features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 6000) 
# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)
filter_genes <- function(seurat_obj) {
  seurat_obj <- subset(seurat_obj, features = integ_features)
  return(seurat_obj)
}
split_seurat_filtered <- lapply(split_seurat, filter_genes)
#perform CCA, find the best buddies or anchors and filter incorrect anchors.
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat_filtered, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
data_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")




