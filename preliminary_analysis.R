#' Authors@R: person("Annice", "Najafi", email = "annicenajafi@tamu.edu")
#' Department of Biomedical Engineering, Texas A&M University, College Station, TX
#' Summer 2024
#' The following program is to perform some preliminary anlaysis on the SHF data

#Load libraries
library(SeuratDisk)
library(SeuratObject)
library(Seurat)
library(tidyverse)
library(AUCell)
library(ggplot2)
library(gridExtra)

#Read relevant genes

read.csv("~/Downloads/WP_HEART_DEVELOPMENT_mouse.csv")->mouse.heart.development
mouse.heart.development$genes

#Read the first file for embryonic day 9.5


#Perform only once to convert to Seurat object/ only E9.5

#UNCOMMENT
#SeuratDisk::Convert("/Users/annicenajafi/Downloads/E9.5_E2S3.MOSTA.h5ad", dest = "h5seurat", overwrite = FALSE)
#Read E9_5 data
E9_5 <- SeuratDisk::LoadH5Seurat("E9.5_E2S3.MOSTA.h5seurat")
#Plot Violin plot to check
Seurat::VlnPlot(E9_5 , features = c("nFeature_count", "nCount_count"), ncol = 2)

#Do QC
E9_5 <- subset(E9_5, nFeature_count > 2000)
E9_5 <- subset(E9_5, nCount_count > 5000)

#Typical procedure
E9_5 <- NormalizeData(E9_5)
E9_5 <- FindVariableFeatures(E9_5)
E9_5 <- ScaleData(E9_5)
E9_5 <- RunPCA(E9_5)


rownames_data <- rownames(E9_5@meta.data)

part_before_underscore <- sub("_.*", "", rownames_data)


# Extract part after underscore using sub
part_after_underscore <- sub(".*_", "", rownames_data)


#my.df <- data.frame(x = part_before_underscore, y = part_after_underscore, expression = FetchData(E9_5, vars = "Tbx5"))


as.matrix(GetAssayData(object = E9_5, slot = "counts"))->E9_5_mat

#c("Tbx5", "Gata4")->genelist
#c("Pitx2", "Tbx1", "Sema3c", "Nkx2.5")->genelist

E9_5 <- FindVariableFeatures(E9_5)

top.5000.hvg <- Seurat::FindVariableFeatures(E9_5, selection.method='vst', nfeatures=5000)

geneSets <-intersect(head(VariableFeatures(top.5000.hvg), 2500), mouse.heart.development$genes)

genelist<-c("Ptch1")
genelist<-c("Isl1", "Mef2c", "Foxh1", "Foxc1", "Foxc2", "Hand2", "Symd1" ,"Fgf8", "Fgf10")
geneSets <- GeneSet(genelist, setName="geneSet1")
# Calculate enrichment scores
AUCell_run(E9_5_mat, geneSets)->AUCell_obj



my.df <- data.frame(x = part_before_underscore, y = part_after_underscore, AUC = t(getAUC(AUCell_obj)), phenotype = E9_5@meta.data$annotation)


ggplot(my.df, aes(x=-as.numeric(y), y=-as.numeric(x), color=geneSet1))+
  geom_point()+scale_color_gradient(low = "#E1AFD1", high = "#B51B75")+
  
  ggtitle("E9.5")+
  labs(colour="Heart\nDevelpment\nGenes")+
  #geom_hline(yintercept =0)+
  theme(
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_blank(),
    #legend.position = "none",
    plot.title = element_text(hjust = 0.5, size=20),
    axis.text = element_blank(),
    text = element_text(size=18))+labs(x="", y="", tag="A")->A


ggplot(my.df, aes(x=-as.numeric(y), y=-as.numeric(x), color=phenotype))+
  geom_point()+
  
  ggtitle("E9.5")+
  labs(colour="Heart\nDevelpment\nGenes")+
  #geom_hline(yintercept =0)+
  theme(
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_blank(),
    #legend.position = "none",
    plot.title = element_text(hjust = 0.5, size=20),
    axis.text = element_blank(),
    text = element_text(size=18))+labs(x="", y="", tag="A")->plt.1



ggplot(my.df, aes(x = -as.numeric(y))) +
  geom_density(color="#B51B75", size=2) +
  theme_classic()+
  labs(x="y")+
  theme(
    plot.title = element_text(hjust = 0.5, size=20),
    text = element_text(size=18))->B


ggplot(my.df, aes(x = -as.numeric(x))) +
  geom_density(color="#B51B75", size=2) +
  theme_classic()+
  coord_flip() +
  labs(x="x")+
  theme(
    plot.title = element_text(hjust = 0.5, size=20),
    text = element_text(size=18))->C



#ggExtra::ggMarginal(A, type = "histogram")

grid.arrange(B, ggplot()+theme_classic(), A, C, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))





#Day E10.5




Raw_data <- Seurat::Read10X(data.dir = '~/Downloads/mice_cardiogenesis_hedgehog/data/E10_5/matrix_files')

read.csv("~/Downloads/mice_cardiogenesis_hedgehog/data/E10_5/metadata.csv")->meta.data

E10_5 <- Seurat::CreateSeuratObject(counts = Raw_data, meta.data = meta.data)



#E10_5[["percent.mt"]] <- Seurat::PercentageFeatureSet(E10_5 , pattern = "^MT-")
Seurat::VlnPlot(E10_5 , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E10_5 <- subset(E10_5, nFeature_RNA > 500)
E10_5 <- subset(E10_5, nCount_RNA > 2000)


E10_5 <- NormalizeData(E10_5)
E10_5 <- FindVariableFeatures(E10_5)
E10_5 <- ScaleData(E10_5)
E10_5 <- RunPCA(E10_5)

features <- c("Gata4", "Gata3")

FeaturePlot(E10_5, features = features)




FetchData(E10_5, vars = "Gata4")

rownames_data <- rownames(E10_5@meta.data)

part_before_underscore <- sub("_.*", "", rownames_data)


# Extract part after underscore using sub
part_after_underscore <- sub(".*_", "", rownames_data)

# 
# my.df <- data.frame(x = part_before_underscore, y = part_after_underscore, expression = FetchData(E10_5, vars = "Tbx5"))
# 
# 
# ggplot(my.df, aes(x=-as.numeric(x), y=-as.numeric(y), color=Tbx5))+
#   geom_point()+scale_color_gradient(low = "#E1AFD1", high = "#B51B75")+
#   ggtitle("E10.5 - TBX5")+
#   #geom_hline(yintercept =0)+
#   theme(
#     # Remove panel border
#     panel.border = element_blank(),
#     # Remove panel grid lines
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     # Remove panel background
#     panel.background = element_blank(),
#     # Add axis line
#     axis.line = element_blank(),
#     #legend.position = "none",
#     plot.title = element_text(hjust = 0.5, size=20),
#     axis.text = element_blank(),
#     text = element_text(size=18))+labs(x="", y="", tag="A")->plt.1
# 
# 
# 
# my.df <- data.frame(x = part_before_underscore, y = part_after_underscore, expression = FetchData(E10_5, vars = "Gata4"))
# 
# 
# ggplot(my.df, aes(x=-as.numeric(x), y=-as.numeric(y), color=Gata4))+
#   geom_point()+scale_color_gradient(low = "#E1AFD1", high = "#B51B75")+
#   ggtitle("E10.5 - GATA4")+
#   #geom_hline(yintercept =0)+
#   theme(
#     # Remove panel border
#     panel.border = element_blank(),
#     # Remove panel grid lines
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     # Remove panel background
#     panel.background = element_blank(),
#     # Add axis line
#     axis.line = element_blank(),
#     #legend.position = "none",
#     plot.title = element_text(hjust = 0.5, size=20),
#     axis.text = element_blank(),
#     text = element_text(size=18))+labs(x="", y="", tag="B")->plt.2
# 
# 
# 
# 
# 
# 
# library(AUCell)
as.matrix(GetAssayData(object = E10_5, slot = "counts"))->E10_5_mat



#c("Tbx5", "Gata4")->genelist
c("Pitx2", "Tbx1", "Sema3c", "Nkx2.5")->genelist
#c("Tbx5", "Gata4")->genelist
c("Pitx2", "Tbx1", "Sema3c", "Nkx2.5", "Gata4", "Tbx5")->genelist

E10_5 <- FindVariableFeatures(E10_5)

top.5000.hvg <- Seurat::FindVariableFeatures(E10_5, selection.method='vst', nfeatures=10000)

geneSets <-intersect(head(VariableFeatures(top.5000.hvg), 6500), mouse.heart.development$genes)

geneSets <- GeneSet(genelist, setName="geneSet1")
# Calculate enrichment scores
AUCell_run(E10_5_mat, geneSets)->AUCell_obj



my.df <- data.frame(x = part_before_underscore, y = part_after_underscore, AUC = t(getAUC(AUCell_obj)), phenotype = E10_5@meta.data$annotation)


ggplot(my.df, aes(x=as.numeric(x), y=-as.numeric(y), color=geneSet1))+
  geom_point()+scale_color_gradient(low = "#E1AFD1", high = "#B51B75")+
  
  ggtitle("E10.5")+
  labs(colour="Heart\nDevelpment\nGenes")+
  #geom_hline(yintercept =0)+
  theme(
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_blank(),
    #legend.position = "none",
    plot.title = element_text(hjust = 0.5, size=20),
    axis.text = element_blank(),
    text = element_text(size=18))+labs(x="", y="", tag="B")->A


ggplot(my.df, aes(x=as.numeric(x), y=-as.numeric(y), color=phenotype))+
  geom_point()+
  
  ggtitle("E10.5")+
  labs(colour="Phenotype")+
  #geom_hline(yintercept =0)+
  theme(
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_blank(),
    #legend.position = "none",
    plot.title = element_text(hjust = 0.5, size=20),
    axis.text = element_blank(),
    text = element_text(size=18))+labs(x="", y="", tag="B")->plt.2


ggplot(my.df, aes(x = as.numeric(x))) +
  geom_density(color="#B51B75", size=2) +
  theme_classic()+
  labs(x="y")+
  theme(
    plot.title = element_text(hjust = 0.5, size=20),
    text = element_text(size=18))->B


ggplot(my.df, aes(x = -as.numeric(y))) +
  geom_density(color="#B51B75", size=2) +
  theme_classic()+
  coord_flip() +
  labs(x="x")+
  theme(
    plot.title = element_text(hjust = 0.5, size=20),
    text = element_text(size=18))->C



#ggExtra::ggMarginal(A, type = "histogram")

grid.arrange(B, ggplot()+theme_classic(), A, C, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))

library(gridExtra)

grid.arrange(B, ggplot(),A, C, 
             widths=c(1, 4), as.table=FALSE, nrow=2)


  #theme(axis.line = element_blank(),
  #      axis.text = element_blank(),
  #      text = element_blank())
  #scale_x_continuous(limits = y_range, expand = expand_scale(mult = 0.1))













#Day E10.5




Raw_data <- Seurat::Read10X(data.dir = '~/Downloads/mice_cardiogenesis_hedgehog/data/E11_5/matrix_files')

read.csv("~/Downloads/mice_cardiogenesis_hedgehog/data/E11_5/metadata.csv")->meta.data

E11_5 <- Seurat::CreateSeuratObject(counts = Raw_data, meta.data = meta.data)



#E10_5[["percent.mt"]] <- Seurat::PercentageFeatureSet(E10_5 , pattern = "^MT-")
Seurat::VlnPlot(E11_5 , features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
E11_5 <- subset(E11_5, nFeature_RNA > 300)
E11_5 <- subset(E11_5, nCount_RNA > 1000)


E11_5 <- NormalizeData(E11_5)
E11_5 <- FindVariableFeatures(E11_5)
E11_5 <- ScaleData(E11_5)
E11_5 <- RunPCA(E11_5)





rownames_data <- rownames(E11_5@meta.data)

part_before_underscore <- sub("_.*", "", rownames_data)

# Extract part after underscore using sub
part_after_underscore <- sub(".*_", "", rownames_data)

# 
# my.df <- data.frame(x = part_before_underscore, y = part_after_underscore, expression = FetchData(E10_5, vars = "Tbx5"))
# 
# 
# ggplot(my.df, aes(x=-as.numeric(x), y=-as.numeric(y), color=Tbx5))+
#   geom_point()+scale_color_gradient(low = "#E1AFD1", high = "#B51B75")+
#   ggtitle("E10.5 - TBX5")+
#   #geom_hline(yintercept =0)+
#   theme(
#     # Remove panel border
#     panel.border = element_blank(),
#     # Remove panel grid lines
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     # Remove panel background
#     panel.background = element_blank(),
#     # Add axis line
#     axis.line = element_blank(),
#     #legend.position = "none",
#     plot.title = element_text(hjust = 0.5, size=20),
#     axis.text = element_blank(),
#     text = element_text(size=18))+labs(x="", y="", tag="A")->plt.1
# 
# 
# 
# my.df <- data.frame(x = part_before_underscore, y = part_after_underscore, expression = FetchData(E10_5, vars = "Gata4"))
# 
# 
# ggplot(my.df, aes(x=-as.numeric(x), y=-as.numeric(y), color=Gata4))+
#   geom_point()+scale_color_gradient(low = "#E1AFD1", high = "#B51B75")+
#   ggtitle("E10.5 - GATA4")+
#   #geom_hline(yintercept =0)+
#   theme(
#     # Remove panel border
#     panel.border = element_blank(),
#     # Remove panel grid lines
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     # Remove panel background
#     panel.background = element_blank(),
#     # Add axis line
#     axis.line = element_blank(),
#     #legend.position = "none",
#     plot.title = element_text(hjust = 0.5, size=20),
#     axis.text = element_blank(),
#     text = element_text(size=18))+labs(x="", y="", tag="B")->plt.2
# 
# 
# 
# 
# 
# 
# library(AUCell)
as.matrix(GetAssayData(object = E11_5, slot = "counts"))->E11_5_mat



#c("Tbx5", "Gata4")->genelist
c("Pitx2", "Tbx1", "Sema3c", "Nkx2.5", "Gata4", "Tbx5")->genelist

E11_5 <- FindVariableFeatures(E11_5)

top.5000.hvg <- Seurat::FindVariableFeatures(E11_5, selection.method='vst', nfeatures=10000)

geneSets <-intersect(head(VariableFeatures(top.5000.hvg), 2500), mouse.heart.development$genes)

geneSets <- GeneSet(genelist, setName="geneSet1")
# Calculate enrichment scores
AUCell_run(E11_5_mat, geneSets)->AUCell_obj



my.df <- data.frame(x = part_before_underscore, y = part_after_underscore, AUC = t(getAUC(AUCell_obj)), phenotype=E11_5@meta.data$annotation)


ggplot(my.df, aes(x=-as.numeric(x), y=-as.numeric(y), color=geneSet1))+
  geom_point()+scale_color_gradient(low = "#E1AFD1", high = "#B51B75")+
  
  ggtitle("E11.5")+
  labs(colour="Heart\nDevelpment\nGenes")+
  #geom_hline(yintercept =0)+
  theme(
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_blank(),
    #legend.position = "none",
    plot.title = element_text(hjust = 0.5, size=20),
    axis.text = element_blank(),
    text = element_text(size=18))+labs(x="", y="", tag="C")->A





ggplot(my.df, aes(x=-as.numeric(x), y=-as.numeric(y), color=phenotype))+
  geom_point()+
  
  ggtitle("E11.5")+
  labs(colour="Phenotype")+
  #geom_hline(yintercept =0)+
  theme(
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_blank(),
    #legend.position = "none",
    plot.title = element_text(hjust = 0.5, size=20),
    axis.text = element_blank(),
    text = element_text(size=18))+labs(x="", y="", tag="C")->plt.3




ggplot(my.df, aes(x = -as.numeric(x))) +
  geom_density(color="#B51B75", size=2) +
  theme_classic()+
  labs(x="y")+
  theme(
    plot.title = element_text(hjust = 0.5, size=20),
    text = element_text(size=18))->B


ggplot(my.df, aes(x = -as.numeric(y))) +
  geom_density(color="#B51B75", size=2) +
  theme_classic()+
  coord_flip() +
  labs(x="x")+
  theme(
    plot.title = element_text(hjust = 0.5, size=20),
    text = element_text(size=18))->C



#ggExtra::ggMarginal(A, type = "histogram")

grid.arrange(B, ggplot()+theme_classic(), A, C, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))




#theme(axis.line = element_blank(),
#      axis.text = element_blank(),
#      text = element_blank())
#scale_x_continuous(limits = y_range, expand = expand_scale(mult = 0.1))



















