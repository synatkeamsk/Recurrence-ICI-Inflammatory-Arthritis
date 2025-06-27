require(Seurat) 
require(tidyverse)
seurat.obj<- readRDS("seurat.obj.final.RDS")

#Dimplot ! 
all.annotated.subsets<- RenameIdents(seurat.obj, `0` ="CD4 Tcm",`11`="Recently activated CD4",  `9`= "PD1hi CXCL13hi CD4", `4`= "Treg", `2` ="Effector CD8", 
                                     `12`= "CD8 Trm", `14`= "Cycing T",`13`= "MAIT",`3` ="Mixed T", `8`= "NK",`7`= "NK T",  `20`= "mDC",  `16`= "CLEC9A+ DCs", 
                                     `17`= " pDCs",  `5`= "Non-inflammatory DCs", `1` ="Classical mono", `10`= "Non-classical mono",`6`= "SPP1+ Mac", `15`="Neutrophil",
                                      `19`= "B cells", `18`= "Syn cells")

## UMAP 
DimPlot(all.annotated.subsets, 
        reduction = "umap", 
        label = FALSE, 
        label.size = 3.5) +
  theme_minimal() + 
  theme(legend.text = element_text(size= 8))
  