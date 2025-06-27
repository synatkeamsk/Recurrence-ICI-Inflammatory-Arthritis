require(tidyverse)
require(pachwork)

#Read the proportion file
proportion<- read.csv("cluster_frequency/cluster_frequency.csv", stringsAsFactors = TRUE)

CD4.ctm<- proportion %>% filter(cluster_name == "CD4 central memory")
CD4.ctm$Group<- factor(CD4.ctm$Group, levels = c("osteoarthritis", "ICI_arthritis"))
CD4.ctm$Arthritis<- factor(CD4.ctm$Arthritis, levels = c("osteoarthritis", "first arthritis", "second arthritis"))

#visualization ! 
# CD4 Tcm
p1<- ggplot(CD4.ctm, aes(x= Group, y= frequency, fill= Group)) +
  theme_classic() +
  geom_boxplot(size= 0.5, alpha= 0.5) +
   geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(color=Arthritis), show.legend = TRUE, size= 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=5.5
  ), axis.title.y = element_blank(), 
  axis.title.x = element_blank(),
  axis.text = element_text(size=8, face= "bold"), 
  axis.text.x = element_blank(), 
  legend.position = "none") + 
  labs(x= "Groups", y= "% of Total Cells", title = "CD4 Tcm")
p1

#Recently activated CD4
ActivatedCD4<- proportion %>% filter(cluster_name == "recently activated CD4")
ActivatedCD4$Group<- factor(ActivatedCD4$Group, levels = c("osteoarthritis", "ICI_arthritis"))
ActivatedCD4$Arthritis<- factor(ActivatedCD4$Arthritis, levels = c("osteoarthritis", "first arthritis", "second arthritis"))

#visualization ! 
p2<- ggplot(ActivatedCD4, aes(x= Group, y= frequency, fill= Group)) +
  theme_classic() +
  geom_boxplot(size= 0.5, alpha= 0.5) +
   geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(color=Arthritis), show.legend = TRUE, size= 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=4.3
  ), axis.title.y = element_blank(), 
  axis.title.x = element_blank(),
  axis.text = element_text(size=8, face = "bold"), 
  axis.text.x = element_blank(), 
  legend.position = "none") + 
  labs(x= "Groups", y= "% Activated CD4 of total cells", title = "recently activated CD4")
p2

# PD1hi CXCL14hi CD4 
TPH<- proportion %>% filter(cluster_name == "PD-1hi CXCL13hi CD4")
TPH$Group<- factor(TPH$Group, levels = c("osteoarthritis", "ICI_arthritis"))
TPH$Arthritis<- factor(TPH$Arthritis, levels = c("osteoarthritis", "first arthritis", "second arthritis"))


#visualization ! 
p3<- ggplot(TPH, aes(x= Group, y= frequency, fill= Group)) +
  theme_classic() +
  geom_boxplot(size= 0.5, alpha= 0.5) +
   geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(color=Arthritis), show.legend = TRUE, size= 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=5.5
  ), axis.title.y = element_blank(), 
  axis.title.x = element_blank(),
  axis.text = element_text(size=8, face = "bold"), 
  axis.text.x = element_blank(),
  legend.position = "none") + 
  labs(x= "Groups", y= "% of Total Cells", title = "PD1hi CXCL13hi CD4")
p3

# MAIT 
MAIT<- proportion %>% filter(cluster_name == "MAIT")
MAIT$Group<- factor(MAIT$Group, levels = c("osteoarthritis", "ICI_arthritis"))
MAIT$Arthritis<- factor(MAIT$Arthritis, levels = c("osteoarthritis", "first arthritis", "second arthritis"))

#visualization ! 
p4<- ggplot(MAIT, aes(x= Group, y= frequency, fill= Group)) +
  theme_classic() +
  geom_boxplot(size= 0.5, alpha= 0.5) +
   geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(color=Arthritis), show.legend = TRUE, size= 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=5.5
  ), axis.title.y = element_blank(), 
  axis.title.x = element_blank(),
  axis.text = element_text(size=8, face = "bold"), 
  axis.text.x = element_blank(), 
  legend.position = "none") + 
  labs(x= "Groups", y= "% MAIT of total cells", title = "MAIT")
p4

# Treg
Treg<- proportion %>% filter(cluster_name == "Treg")
Treg$Group<- factor(Treg$Group, levels = c("osteoarthritis", "ICI_arthritis"))
Treg$Arthritis<- factor(Treg$Arthritis, levels = c("osteoarthritis", "first arthritis", "second arthritis"))


#visualization ! 
p5<- ggplot(Treg, aes(x= Group, y= frequency, fill= Group)) +
  theme_classic() +
  geom_boxplot(size= 0.5, alpha= 0.5) +
   geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(color=Arthritis), show.legend = TRUE, size= 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=5.5
  ), axis.title.y = element_blank(), 
  axis.title.x = element_blank(),
  axis.text = element_text(size=8, face= "bold"), 
  axis.text.x = element_blank(),
  legend.position = "none") + 
  labs(x= "Groups", y= "% Treg of total cells", title = "Treg")
p5

# Effector CD8 
effector.cd8<- proportion %>% filter(cluster_name == "Effector CD8")
effector.cd8$Group<- factor(effector.cd8$Group, levels = c("osteoarthritis", "ICI_arthritis"))
effector.cd8$Arthritis<- factor(effector.cd8$Arthritis, levels = c("osteoarthritis", "first arthritis", "second arthritis"))


#visualization ! 
p6<- ggplot(effector.cd8, aes(x= Group, y= frequency, fill= Group)) +
  theme_classic() +
  geom_boxplot(size= 0.5, alpha= 0.5) +
   geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(color=Arthritis), show.legend = TRUE, size= 2.5) +
 theme(plot.title = element_text(hjust = 0.5, face = "bold", size=5.5
  ), axis.title.y = element_blank(), 
  axis.title.x = element_blank(),
  axis.text = element_text(size=8, face = "bold"), 
  axis.text.x = element_blank(),
  legend.position = "none") + 
  labs(x= "Groups", y= "% of Total Cells", title = "Effector CD8")
p6

# CD8 Trm
TRM.cd8<- proportion %>% filter(cluster_name == "Trm effector CD8 ")
TRM.cd8$Group<- factor(TRM.cd8$Group, levels = c("osteoarthritis", "ICI_arthritis"))
TRM.cd8$Arthritis<- factor(TRM.cd8$Arthritis, levels = c("osteoarthritis", "first arthritis", "second arthritis"))


#visualization ! 
p7<- ggplot(TRM.cd8, aes(x= Group, y= frequency, fill= Group)) +
  theme_classic() +
  geom_boxplot(size= 0.5, alpha= 0.5) +
   geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(color=Arthritis), show.legend = TRUE, size= 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=5.5
  ), axis.title.y = element_blank(), 
  axis.title.x = element_blank(),
  axis.text = element_text(size=8, face= "bold"), 
  axis.text.x = element_blank(),
  legend.position = "none") + 
  labs(x= "Groups", y= "% Trm effector CD8 of total cells", title = "Trm CD8")
p7

# Cycling T cell
cyclingT<- proportion %>% filter(cluster_name == "Cycing T")
cyclingT$Group<- factor(cyclingT$Group, levels = c("osteoarthritis", "ICI_arthritis"))
cyclingT$Arthritis<- factor(cyclingT$Arthritis, levels = c("osteoarthritis", "first arthritis", "second arthritis"))


#visualization ! 
p8<- ggplot(cyclingT, aes(x= Group, y= frequency, fill= Group)) +
  theme_classic() +
  geom_boxplot(size= 0.5, alpha= 0.5) +
   geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(color=Arthritis), show.legend = TRUE, size= 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=5.5
  ), axis.title.y = element_blank(), 
  axis.title.x = element_blank(),
  axis.text = element_text(size=8, face= "bold"), 
  axis.text.x = element_blank(),
  legend.position = "none") + 
  labs(x= "Groups", y= "% cycling T of total cells", title = "Cycling T")
p8

## NK T cell
NKT<- proportion %>% filter(cluster_name == "NKT")
NKT$Group<- factor(NKT$Group, levels = c("osteoarthritis", "ICI_arthritis"))
NKT$Arthritis<- factor(NKT$Arthritis, levels = c("osteoarthritis", "first arthritis", "second arthritis"))


#visualization ! 
p9<- ggplot(NKT, aes(x= Group, y= frequency, fill= Group)) +
   theme_classic() +
  geom_boxplot(size= 0.5, alpha= 0.5) +
   geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(color=Arthritis), show.legend = TRUE, size= 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=5.5
  ), axis.title.y = element_blank(), 
  axis.title.x = element_blank(),
  axis.text = element_text(size=8, face= "bold"), 
  axis.text.x = element_blank(), 
  legend.position = "none") + 
  labs(x= "Groups", y= "% NK T of total cells", title = "NK T")
p9

# NK 
Nk<- proportion %>% filter(cluster_name == "NK")
Nk$Group<- factor(Nk$Group, levels = c("osteoarthritis", "ICI_arthritis"))
Nk$Arthritis<- factor(Nk$Arthritis, levels = c("osteoarthritis", "first arthritis", "second arthritis"))

#visualization ! 
p11<- ggplot(Nk, aes(x= Group, y= frequency, fill= Group)) +
  theme_classic() +
  geom_boxplot(size= 0.5, alpha= 0.5) +
   geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(color=Arthritis), show.legend = TRUE, size= 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=5.5
  ), axis.title.y = element_blank(), 
  axis.title.x = element_blank(),
  axis.text = element_text(size=8, face= "bold"), 
  axis.text.x = element_blank(), 
  legend.position = "none") + 
  labs(x= "Groups", y= "% of Totol Cells", title = "NK")
p11

# B cell
Bcell<- proportion %>% filter(cluster_name == "B cells")
Bcell$Arthritis<- factor(Bcell$Arthritis, levels = c("osteoarthritis","first arthritis", "second arthritis"))
Bcell$Group<- factor(Bcell$Group, levels = c("osteoarthritis", "ICI_arthritis"))


#visualization ! 
p12<- ggplot(Bcell, aes(x= Group, y= frequency, fill= Group)) +
  theme_classic() +
  geom_boxplot(size= 0.5, alpha= 0.5) +
   geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(color=Arthritis), show.legend = TRUE, size= 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=5.5
  ), axis.title.y = element_blank(), 
  axis.title.x = element_blank(),
  axis.text = element_text(size=8, face= "bold"), 
  axis.text.x = element_blank(), 
  legend.position = "none") +
  labs(x= "Groups", y= "% B cell of total cells", title = "B cell")
p12

# Classical monocyte
classmono<- proportion %>% filter(cluster_name == "Classical Monocytes")
classmono$Arthritis<- factor(classmono$Arthritis, levels = c("osteoarthritis","first arthritis", "second arthritis"))
classmono$Group<- factor(classmono$Group, levels = c("osteoarthritis", "ICI_arthritis"))

#visualization ! 
p13<- ggplot(mono.macro, aes(x= Group, y= frequency, fill= Group)) +
  theme_classic() +
  geom_boxplot(size= 0.5, alpha= 0.5) +
   geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(color=Arthritis), show.legend = TRUE, size= 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=5.5
  ), axis.title.y = element_blank(), 
  axis.title.x = element_blank(),
  axis.text = element_text(size=8, face = "bold"), 
  axis.text.x = element_blank(), 
  legend.position = "none") + 
  labs(x= "Groups", y= "% classical mono of total cells", title = "Classical mono")
p13

#Neutrophil
neutro<- proportion %>% filter(cluster_name == "Non-classical monocytes")
neutro$Arthritis<- factor(neutro$Arthritis, levels = c("osteoarthritis","first arthritis", "second arthritis"))
neutro$Group<- factor(neutro$Group, levels = c("osteoarthritis", "ICI_arthritis"))


#visualization ! 
p14<- ggplot(neutro, aes(x= Group, y= frequency, fill= Group)) +
  theme_classic() +
  geom_boxplot(size= 0.5, alpha= 0.5) +
   geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(color=Arthritis), show.legend = TRUE, size= 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=5.4
  ), axis.title.y = element_blank(), 
  axis.title.x = element_blank(),
  axis.text = element_text(size=8, face= "bold"), 
  axis.text.x = element_blank(), 
  legend.position = "none") + 
  labs(x= "Groups", y= "% Non-class mono of all cells", title = "Non-classical mono")
p14

#SPP1+ Macrophages
SPP1macrophage<- proportion %>% filter(cluster_name == "SPP1+ Macrophages")
SPP1macrophage$Arthritis<- factor(SPP1macrophage$Arthritis, levels = c("osteoarthritis","first arthritis", "second arthritis"))
SPP1macrophage$Group<- factor(SPP1macrophage$Group, levels = c("osteoarthritis", "ICI_arthritis"))


#visualization ! 
p15<- ggplot(SPP1macrophage, aes(x= Group, y= frequency, fill= Group)) +
  theme_classic() +
  geom_boxplot(size= 0.5, alpha= 0.5) +
   geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(color=Arthritis), show.legend = TRUE, size= 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=5.5
  ), axis.title.y = element_blank(), 
  axis.title.x = element_blank(),
  axis.text = element_text(size=8, face= "bold"), 
  axis.text.x = element_blank(), 
  legend.position = "none") + 
  labs(x= "Groups", y= "% SPP1+ Macrophages of CD45", title = "SPP1+ Macrophages")
p15

#CD1c+ Non-inflammatory DCs
monocyte<- proportion %>% filter(cluster_name == "CD1c+ Non-inflammatory DCs")
monocyte$Arthritis<- factor(monocyte$Arthritis, levels = c("osteoarthritis","first arthritis", "second arthritis"))
monocyte$Group<- factor(monocyte$Group, levels = c("osteoarthritis", "ICI_arthritis"))

#visualization ! 
p16<- ggplot(monocyte, aes(x= Group, y= frequency, fill= Group)) +
  theme_classic() +
  geom_boxplot(size= 0.5, alpha= 0.5) +
   geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(color=Arthritis), show.legend = TRUE, size= 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=5.4
  ), axis.title.y = element_blank(), 
  axis.title.x = element_blank(),
  axis.text = element_text(size=8, face= "bold"), 
  axis.text.x = element_blank(), 
  legend.position = "none") + 
  labs(x= "Groups", y= "% of Total Cells", title = "Non-inflamed DC")
p16

# CLEC9A+ DCs
CLEC9ADC<- proportion %>% filter(cluster_name == "CLEC9A+ DCs")
CLEC9ADC$Arthritis<- factor(CLEC9ADC$Arthritis, levels = c("osteoarthritis","first arthritis", "second arthritis"))
CLEC9ADC$Group<- factor(CLEC9ADC$Group, levels = c("osteoarthritis", "ICI_arthritis"))


#visualization ! 
p17<- ggplot(CLEC9ADC, aes(x= Group, y= frequency, fill= Group)) +
  theme_classic() +
  geom_boxplot(size= 0.5, alpha= 0.5) +
   geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(color=Arthritis), show.legend = TRUE, size= 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=5.5
  ), axis.title.y = element_blank(), 
  axis.title.x = element_blank(),
  axis.text = element_text(size=8, face = "bold"), 
  axis.text.x = element_blank(), 
  legend.position = "none") + 
  labs(x= "Groups", y= "% CLEC9A+ DCs of all cells", title = "CLEC9A+ DCs")
p17

#pDC
pDC<- proportion %>% filter(cluster_name == "pDC")
pDC$Arthritis<- factor(pDC$Arthritis, levels = c("osteoarthritis","first arthritis", "second arthritis"))
pDC$Group<- factor(pDC$Group, levels = c("osteoarthritis", "ICI_arthritis"))


#visualization ! 
p18<- ggplot(pDC, aes(x= Group, y= frequency, fill= Group)) +
  theme_classic() +
  geom_boxplot(size= 0.5, alpha= 0.5) +
   geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(color=Arthritis), show.legend = TRUE, size= 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=5.5
  ), axis.title.y = element_blank(), 
  axis.title.x = element_blank(),
  axis.text = element_text(size=8, face= "bold"), 
  axis.text.x = element_blank(), 
  legend.position = "none") + 
  labs(x= "Groups", y= "% of Total Cells", title = "pDC")
p18

#mDC 
mdC<- proportion %>% filter(cluster_name == "mDC")
mdC$Arthritis<- factor(mdC$Arthritis, levels = c("osteoarthritis","first arthritis", "second arthritis"))
mdC$Group<- factor(mdC$Group, levels = c("osteoarthritis", "ICI_arthritis"))

#visualization ! 
p19<- ggplot(mdC, aes(x= Group, y= frequency, fill= Group)) +
  theme_classic() +
  geom_boxplot(size= 0.5, alpha= 0.5) +
   geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(color=Arthritis), show.legend = TRUE, size= 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=5.5
  ), axis.title.y = element_blank(), 
  axis.title.x = element_blank(),
  axis.text = element_text(size=8, face= "bold"), 
  axis.text.x = element_blank(), 
  legend.position = "none") + 
  labs(x= "Groups", y= "% mDC of CD45", title = "mDC")
p19

#Neutrophil 
Neutrophil<- proportion %>% filter(cluster_name == "Neutrophil")
Neutrophil$Arthritis<- factor(Neutrophil$Arthritis, levels = c("osteoarthritis","first arthritis", "second arthritis"))
Neutrophil$Group<- factor(Neutrophil$Group, levels = c("osteoarthritis", "ICI_arthritis"))

#visualization ! 
p20<- ggplot(Neutrophil, aes(x= Group, y= frequency, fill= Group)) +
  theme_classic() +
  geom_boxplot(size= 0.5, alpha= 0.5) +
   geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(color=Arthritis), show.legend = TRUE, size= 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=5.5
  ), axis.title.y = element_blank(), 
  axis.title.x = element_blank(),
  axis.text = element_text(size=8, face= "bold"), 
  axis.text.x = element_blank(), 
  legend.position = "none") + 
  labs(x= "Groups", y= "% Total cells", title = "Neutrophil")
p20

#Synovial cells
synovial<- proportion %>% filter(cluster_name == "Synovial cells")
synovial$Arthritis<- factor(synovial$Arthritis, levels = c("osteoarthritis","first arthritis", "second arthritis"))
synovial$Group<- factor(synovial$Group, levels = c("osteoarthritis", "ICI_arthritis"))

#visualization ! 
p21<- ggplot(synovial, aes(x= Group, y= frequency, fill= Group)) +
  theme_classic() +
  geom_boxplot(size= 0.5, alpha= 0.5) +
   geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(color=Arthritis), show.legend = TRUE, size= 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=5.5
  ), axis.title.y = element_blank(), 
  axis.title.x = element_blank(),
  axis.text = element_text(size=8, face = "bold"), 
  axis.text.x = element_blank(), 
  legend.position = "none", 
  legend.text = element_text(face = "bold"), 
  legend.title = element_text(face = "bold")) + 
  labs(x= "Groups", y= "% Synovial cells of all cells", title = "Synovial cells")
p21

# combine each plots 
proportion1<- (p20|p13|p14|p15|p19|p18|p16|p17|p1|p2)
ggsave("immune_fraction.new.pdf", plot = proportion1, height=1.5, width = 9, units = "in", dpi = 300)

proportion2<- (p3|p5|p6|p7|p8|p4|p11|p9|p12|p21)
ggsave("immune_fraction1.new.pdf", plot = proportion2, height=1.5, width = 9, units = "in", dpi = 300)








































































