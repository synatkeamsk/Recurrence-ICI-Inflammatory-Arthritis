require(Seurat) 
require(tidyverse)
require(viridis)

#Read the annotated obj
all.annotated.subsets<- readRDS("all.clusters.annotated.obj.rds")

#Gene list !! 
genes <- c("CD3D", "MKI67", "CD4", "IL2RA", "FOXP3", "CD8A", "TRDV2", "IL7R", "KLRG1", "NCAM1", "CD19", "TBX21", "MS4A1", "CD27", "CD38", "SDC1", "CD14", 
"FCGR3A", "S100A8", "S100A9", "HLA-DRB1", "CD1C", "CXCR3", "CXCR6", "CCR6", "CCR7", "CXCL13", "TOX", "LAG3", "PDCD1","CTLA4", "HAVCR2", "TGIT", "EOMES", 
"GZMB", "GZMH", "GZMK", "IFNG", "TNF", "IL10", "IL17A", "IL21", "TGFB1", "IL2")

# Dotplot 
DotPlot(all.annotated.subsets, features = genes, dot.scale = 4.5) + 
 geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="inferno") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15), 
        axis.title.y= element_blank(), 
        axis.text = element_text(size=10, face = "bold"), 
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.2, face = "bold", size= 10), 
        axis.title.x = element_blank(),
        legend.title = element_text(size= 12, face = "bold"), 
        legend.text = element_text(size= 10, face = "bold"), 
        legend.position = "right") 