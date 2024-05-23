#' Authors@R: person("Annice", "Najafi", email = "annicenajafi@tamu.edu")
#' Department of Biomedical Engineering, Texas A&M University, College Station, TX
#' Summer 2024
#' The following program is to separate out the SHF to OFT trajectory and extract
#' the expression of ligand and receptor genes related to the hedgehog pathway



#Load libraries
library(SeuratDisk)
library(SeuratObject)
library(Seurat)
library(tidyverse)
library(AUCell)
library(ggplot2)
library(gridExtra)


#Read E10.5 data

Raw_data <- Seurat::Read10X(data.dir = '~/Downloads/mice_cardiogenesis_hedgehog/data/E10_5/matrix_files')

read.csv("~/Downloads/mice_cardiogenesis_hedgehog/data/E10_5/metadata.csv")->meta.data

E10_5 <- Seurat::CreateSeuratObject(counts = Raw_data, meta.data = meta.data)


#Follow the typical steps for QC

#E10_5[["percent.mt"]] <- Seurat::PercentageFeatureSet(E10_5 , pattern = "^MT-")
Seurat::VlnPlot(E10_5 , features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
E10_5 <- subset(E10_5, nFeature_RNA > 500)
E10_5 <- subset(E10_5, nCount_RNA > 2000)


E10_5 <- NormalizeData(E10_5)
E10_5 <- FindVariableFeatures(E10_5)
E10_5 <- ScaleData(E10_5)
E10_5 <- RunPCA(E10_5)


#Extract the locations
rownames_data <- rownames(E10_5@meta.data)
part_before_underscore <- sub("_.*", "", rownames_data)
# Extract part after underscore using sub
part_after_underscore <- sub(".*_", "", rownames_data)


#convert to matrix to use later
as.matrix(GetAssayData(object = E10_5, slot = "counts"))->E10_5_mat


#Extract the trajectory based on the expression of these two genes
c("Ptch1", "Osr1")->genelist

geneSets <- GeneSet(genelist, setName="geneSet1")
# Calculate enrichment scores
AUCell_run(E10_5_mat, geneSets)->AUCell_obj
enrichment.df <- t(getAUC(AUCell_obj))


#Initial plot to visualize the enrichment of the two biomarkers distribution
ggplot(enrichment.df, aes(x = geneSet1))+geom_density(size=2, color="#AD88C6")+
  ggtitle("E10.5 - Osr1/Ptch1 Enrichment") +
  labs(colour="Ptch1") +
  geom_vline(xintercept = 0.04, linetype="dotted", size=1.5, color="#E1AFD1")+
  theme_classic()+
  theme(
    #panel.border = element_blank(),
    #panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    #axis.line = element_blank(),
    plot.title = element_text(hjust = 0.5, size=20),
    #axis.text = element_blank(),
    text = element_text(size=18)
  ) + 
  labs(x="", y="", tag="A")->plt.density.A


#Use a cutoff of 0.4 to isolate cells that express Osr1 and Ptch1
E10_5@meta.data$OP_enrichment <- enrichment.df[,1]

E10_5@meta.data$OP_enrichment[E10_5@meta.data$OP_enrichment > 0.04] <- "Y"
E10_5@meta.data$OP_enrichment[E10_5@meta.data$OP_enrichment <= 0.04] <- "N"

#Extract the locations
my.df <- data.frame(x = part_before_underscore, y = part_after_underscore, OP_enrich = E10_5@meta.data$OP_enrichment,  AUC = t(getAUC(AUCell_obj)), phenotype = E10_5@meta.data$annotation)

  


#We are only interested in the cells near the heart
E10_5@meta.data$x_loc <- part_before_underscore
E10_5@meta.data$y_loc <- part_after_underscore


#Get the center of the heart 
circle_center_x_loc <- mean(as.numeric(E10_5@meta.data$x_loc[E10_5@meta.data$annotation == "Heart"]))
circle_center_y_loc <- mean(-as.numeric(E10_5@meta.data$y_loc[E10_5@meta.data$annotation == "Heart"]))


#Extract only regions enclosed within a radius of 30 to 45 from the center of the circle
E10_5@meta.data$refined_OP_path<-"N"
E10_5@meta.data[sqrt((as.numeric(E10_5@meta.data$x_loc) - circle_center_x_loc)^2 + (-as.numeric(E10_5@meta.data$y_loc) - circle_center_y_loc)^2)<45 & sqrt((as.numeric(E10_5@meta.data$x_loc) - circle_center_x_loc)^2 + (-as.numeric(E10_5@meta.data$y_loc) - circle_center_y_loc)^2)>30,]$refined_OP_path<-"Y"



my.df <- data.frame(x = part_before_underscore, y = part_after_underscore, OP_enrich = E10_5@meta.data$OP_enrichment, OP_path = E10_5@meta.data$refined_OP_path,  AUC = t(getAUC(AUCell_obj)), phenotype = E10_5@meta.data$annotation)


ggplot(my.df, aes(x=as.numeric(x), y=-as.numeric(y), color=OP_path))+
  geom_point()+
  #scale_color_gradient(low = "#E1AFD1", high = "#B51B75")+
  scale_color_manual(values = c("#E1AFD1", "#AD88C6"))+
  ggtitle("E10.5")+
  labs(colour="Heart vicinity")+
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
    text = element_text(size=18))+labs(x="", y="", tag="B")+
  geom_point(aes(x = mean(as.numeric(E10_5@meta.data$x_loc[E10_5@meta.data$annotation == "Heart"])), y = -mean(as.numeric(E10_5@meta.data$y_loc[E10_5@meta.data$annotation == "Heart"]))), color = "blue", size = 2) +
  annotate("path",
           x=mean(as.numeric(E10_5@meta.data$x_loc[E10_5@meta.data$annotation == "Heart"]))+45*cos(seq(0,2*3.14,length.out=100)),
           y=-mean(as.numeric(E10_5@meta.data$y_loc[E10_5@meta.data$annotation == "Heart"]))+45*sin(seq(0,2*3.14,length.out=100)))+
  annotate("path",
           x=mean(as.numeric(E10_5@meta.data$x_loc[E10_5@meta.data$annotation == "Heart"]))+30*cos(seq(0,2*3.14,length.out=100)),
           y=-mean(as.numeric(E10_5@meta.data$y_loc[E10_5@meta.data$annotation == "Heart"]))+30*sin(seq(0,2*3.14,length.out=100)))->plt.heart.vicinity



#

#function to find the Euclidean distance
euclidean_distance <- function(x1, y1, x2, y2) {
  sqrt((x1 - x2)^2 + (y1 - y2)^2)
}


as.numeric(E10_5@meta.data$x_loc)->E10_5@meta.data$x_loc
as.numeric(E10_5@meta.data$y_loc)->E10_5@meta.data$y_loc
# Initialize a column to store the average number of "Y" values for each cell
E10_5@meta.data$avg_Y_neighbors <- 0

# Loop through each cell to find its top 10 nearest neighbors and calculate the average number of "Y" values
for (i in 1:nrow(E10_5@meta.data)) {
  # Calculate distances from the current cell to all other cells
  distances <- mapply(euclidean_distance, E10_5@meta.data$x_loc[i], -E10_5@meta.data$y_loc[i], E10_5@meta.data$x_loc, -E10_5@meta.data$y_loc)
  
  # Exclude the distance to itself by setting it to a high value
  distances[i] <- 100000
  
  # Find the indices of the top 10 nearest neighbors
  neighbor_indices <- order(distances)[1:10]
  
  # Calculate the average number of "Y" values among the nearest neighbors
  avg_Y <- mean(E10_5@meta.data$OP_enrich[neighbor_indices] == "Y")
  
  # Store the average in the dataframe
  E10_5@meta.data$avg_Y_neighbors[i] <- avg_Y
}

my.df$y_neigh <- E10_5@meta.data$avg_Y_neighbors



#my.df$heart<-"N"
#my.df[my.df$phenotype=="Heart",]->my.df$heart<-"Y"

#For demonstration purposes plot the enrichment path
ggplot()+
  geom_point(data = my.df, aes(x=as.numeric(x), y=-as.numeric(y), color=y_neigh, stroke=0.01))+
  #scale_color_gradient(low = "#E1AFD1", high = "#B51B75")+
  #scale_color_manual(values = c("#E1AFD1", "#AD88C6"))+
  scale_shape_manual(values = c(19, 8))+
  ggtitle("E10.5")+
  labs(colour="Enriched Neighbors")+
  #geom_hline(yintercept =0)+
  theme_classic()+
  theme(
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    #legend.position = "none",
    plot.title = element_text(hjust = 0.5, size=20),
    axis.text = element_blank(),
    text = element_text(size=18))+labs(x="", y="", tag="A")+
  geom_point(aes(x = mean(as.numeric(E10_5@meta.data$x_loc[E10_5@meta.data$annotation == "Heart"])), y = -mean(as.numeric(E10_5@meta.data$y_loc[E10_5@meta.data$annotation == "Heart"]))), color = "darkred", size = 2) +
  annotate("path",
           x=mean(as.numeric(E10_5@meta.data$x_loc[E10_5@meta.data$annotation == "Heart"]))+45*cos(seq(0,2*3.14,length.out=100)),
           y=-mean(as.numeric(E10_5@meta.data$y_loc[E10_5@meta.data$annotation == "Heart"]))+45*sin(seq(0,2*3.14,length.out=100)))+
  annotate("path",
           x=mean(as.numeric(E10_5@meta.data$x_loc[E10_5@meta.data$annotation == "Heart"]))+30*cos(seq(0,2*3.14,length.out=100)),
           y=-mean(as.numeric(E10_5@meta.data$y_loc[E10_5@meta.data$annotation == "Heart"]))+30*sin(seq(0,2*3.14,length.out=100)))





#In this case we are inverting y coordinates so would have to take max value for min loc and vice versa

min_locs_on_path <- E10_5@meta.data[E10_5@meta.data$y_loc == max(E10_5@meta.data[E10_5@meta.data$avg_Y_neighbors>0.5 & E10_5@meta.data$refined_OP_path=="Y",]$y_loc) & E10_5@meta.data$refined_OP_path == "Y" & E10_5@meta.data$final_traj == "Y",][1,]
max_locs_on_path <- E10_5@meta.data[E10_5@meta.data$y_loc == min(E10_5@meta.data[E10_5@meta.data$avg_Y_neighbors>0.5 & E10_5@meta.data$refined_OP_path=="Y",]$y_loc) & E10_5@meta.data$refined_OP_path == "Y" & E10_5@meta.data$final_traj == "Y",][1,]




E10_5@meta.data$final_traj<-"N" #Ensure it is on the trajectory
#E10_5@meta.data[y_neigh>0.5]->E10_5@meta.data$final_traj
E10_5@meta.data[E10_5@meta.data$avg_Y_neighbors>0.5,]$final_traj<-"Y" #Filter for neighbor enrichment

E10_5@meta.data$final_traj -> my.df$final_traj



ggplot()+
  geom_point(data = my.df, aes(x=as.numeric(x), y=-as.numeric(y), color=final_traj, stroke=0.01))+
  #scale_color_gradient(low = "#E1AFD1", high = "#B51B75")+
  #scale_color_manual(values = c("#E1AFD1", "#AD88C6"))+
  scale_shape_manual(values = c(19, 8))+
  ggtitle("E10.5")+
  labs(colour="Enriched Neighbors")+
  #geom_hline(yintercept =0)+
  theme_classic()+
  theme(
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    #legend.position = "none",
    plot.title = element_text(hjust = 0.5, size=20),
    axis.text = element_blank(),
    text = element_text(size=18))+labs(x="", y="", tag="A")+
  geom_point(aes(x = mean(as.numeric(E10_5@meta.data$x_loc[E10_5@meta.data$annotation == "Heart"])), y = -mean(as.numeric(E10_5@meta.data$y_loc[E10_5@meta.data$annotation == "Heart"]))), color = "darkred", size = 2) +
  annotate("path",
           x=mean(as.numeric(E10_5@meta.data$x_loc[E10_5@meta.data$annotation == "Heart"]))+45*cos(seq(0,2*3.14,length.out=100)),
           y=-mean(as.numeric(E10_5@meta.data$y_loc[E10_5@meta.data$annotation == "Heart"]))+45*sin(seq(0,2*3.14,length.out=100)))+
  annotate("path",
           x=mean(as.numeric(E10_5@meta.data$x_loc[E10_5@meta.data$annotation == "Heart"]))+30*cos(seq(0,2*3.14,length.out=100)),
           y=-mean(as.numeric(E10_5@meta.data$y_loc[E10_5@meta.data$annotation == "Heart"]))+30*sin(seq(0,2*3.14,length.out=100)))



E10_5@meta.data[(-E10_5@meta.data$y_loc> -min_locs_on_path$y_loc) & (-max_locs_on_path$y_loc > -E10_5@meta.data$y_loc) & E10_5@meta.data$refined_OP_path == "Y" ,]$final_traj <- "Y"
E10_5@meta.data[(E10_5@meta.data$x_loc < max_locs_on_path$x_loc) & E10_5@meta.data$final_traj=="Y" & 
                  E10_5@meta.data$y_loc> max_locs_on_path$y_loc + (-max_locs_on_path$y_loc - circle_center_y_loc - 30), ]$final_traj<-"N"

E10_5@meta.data[E10_5@meta.data$refined_OP_path=="N",]$final_traj<-"N" #ensure it is within the two circles


my.df <- data.frame(x = part_before_underscore, y = part_after_underscore, OP_enrich = E10_5@meta.data$OP_enrichment, OP_path = E10_5@meta.data$refined_OP_path,  AUC = t(getAUC(AUCell_obj)), phenotype = E10_5@meta.data$annotation, y_neigh = E10_5@meta.data$avg_Y_neighbors, final_traj = E10_5@meta.data$final_traj)


ggplot()+
  geom_point(data = my.df, aes(x=as.numeric(x), y=-as.numeric(y), color=final_traj, stroke=0.01))+
  #scale_color_gradient(low = "#E1AFD1", high = "#B51B75")+
  #scale_color_manual(values = c("#E1AFD1", "#AD88C6"))+
  scale_shape_manual(values = c(19, 8))+
  ggtitle("E10.5")+
  labs(colour="Enriched Neighbors")+
  #geom_hline(yintercept =0)+
  theme_classic()+
  theme(
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    #legend.position = "none",
    plot.title = element_text(hjust = 0.5, size=20),
    axis.text = element_blank(),
    text = element_text(size=18))+labs(x="", y="", tag="A")+
  geom_point(aes(x = mean(as.numeric(E10_5@meta.data$x_loc[E10_5@meta.data$annotation == "Heart"])), y = -mean(as.numeric(E10_5@meta.data$y_loc[E10_5@meta.data$annotation == "Heart"]))), color = "darkred", size = 2) +
  geom_point(aes(x= min_locs_on_path$x_loc, y=-min_locs_on_path$y_loc), color="yellow")+
  geom_point(aes(x= max_locs_on_path$x_loc, y=-max_locs_on_path$y_loc), color="yellow")+
  annotate("path",
           x=mean(as.numeric(E10_5@meta.data$x_loc[E10_5@meta.data$annotation == "Heart"]))+45*cos(seq(0,2*3.14,length.out=100)),
           y=-mean(as.numeric(E10_5@meta.data$y_loc[E10_5@meta.data$annotation == "Heart"]))+45*sin(seq(0,2*3.14,length.out=100)))+
  annotate("path",
           x=mean(as.numeric(E10_5@meta.data$x_loc[E10_5@meta.data$annotation == "Heart"]))+30*cos(seq(0,2*3.14,length.out=100)),
           y=-mean(as.numeric(E10_5@meta.data$y_loc[E10_5@meta.data$annotation == "Heart"]))+30*sin(seq(0,2*3.14,length.out=100)))




#Center and radii of the circles
center_x <- circle_center_x_loc
center_y <- circle_center_y_loc
inner_radius <- 30
outer_radius <- 45




marker_expression_SHF_OFT_traj<- function(sliding_window_size = 20, center_x, center_y, inner_radius, outer_radius, biomarker){
  
  starting_angle <- - acos(((min_locs_on_path$x_loc-center_x)/inner_radius)) # starting SHF
  
  #ending_angle <- asin(-(max_locs_on_path$y_loc+center_y)/outer_radius) #ending OFT
  ending_angle <- acos((max_locs_on_path$x_loc-center_x)/inner_radius)
  
  sliding_x <- c()
  sliding_y <- c()
  
  for(k in 1:sliding_window_size){
    # Define the angle in radians
    midpoint_angle <- starting_angle + (ending_angle-starting_angle)/sliding_window_size*k
    
    # Calculate the midpoint coordinates on the inner and outer circles
    x_midpoint <- center_x + ((outer_radius-inner_radius)/2 + inner_radius) * cos(midpoint_angle)
    y_midpoint <- center_y + ((outer_radius-inner_radius)/2 + inner_radius)* sin(midpoint_angle)
    
    sliding_x[k]<-x_midpoint
    sliding_y[k]<-y_midpoint
    
  }
  
  windows <- data.frame(x=sliding_x, y=sliding_y)
  
  ligand.df.orig <- FetchData(E10_5, vars = biomarker) 
  sliding_ligand.df <- list()
  
  for(k in 1:sliding_window_size){
    
    max_x <- sliding_x[k] + (outer_radius-inner_radius)/2
    max_y <- sliding_y[k] + (outer_radius-inner_radius)/2
    min_x <- sliding_x[k] - (outer_radius-inner_radius)/2
    min_y <- sliding_y[k] - (outer_radius-inner_radius)/2
    
    ligand.df.orig->ligand.df
    colnames(ligand.df)<-c("ligand")
    # Parse rownames to extract x and y coordinates
    ligand.df$x <- as.numeric(sapply(strsplit(rownames(ligand.df), "_"), `[`, 1))
    ligand.df$y <- as.numeric(sapply(strsplit(rownames(ligand.df), "_"), `[`, 2))
    
    
    ligand.df[(ligand.df$x<max_x & ligand.df$x>min_x & -ligand.df$y>min_y & -ligand.df$y<max_y),]->sliding_ligand.df[[k]]
    
  }
  
  return(list(sliding_ligand.df, windows))
}



sliding_ligand.df <- marker_expression_SHF_OFT_traj(sliding_window_size = 20, center_x, center_y, inner_radius, outer_radius, "Ptch1")

sliding_ligand.df[[2]]->windows
colnames(windows)<-c("sliding_x", "sliding_y")
sliding_ligand.df[[1]]->sliding_ligand.df

ggplot()+
  geom_point(data=my.df, aes(x=as.numeric(x), y=-as.numeric(y), color=final_traj))+
  #geom_rect(data = slices_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
  #          fill = "red", alpha = 0.9, inherit.aes = F)+
  #scale_color_gradient(low = "#E1AFD1", high = "#B51B75")+
  #scale_color_manual(values = c("#E1AFD1", "#AD88C6"))+
  ggtitle("Trajectory")+
  labs(colour="Final\nTrajectory")+
  #geom_vline(xintercept=max_locs_on_path$x)+
  #geom_hline(yintercept=-max_locs_on_path$y)+
  #geom_vline(xintercept=min_locs_on_path$x)+
  #geom_hline(yintercept=-min_locs_on_path$y)+
  scale_color_manual(values = c("#E1AFD1", "#AD88C6"))+
  theme_classic()+
  #geom_hline(yintercept =0)+
  theme(
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    #legend.position = "none",
    plot.title = element_text(hjust = 0.5, size=20),
    axis.text = element_blank(),
    text = element_text(size=18))+labs(x="", y="", tag="B")+
  geom_point(aes(x= min_locs_on_path$x_loc, y=-min_locs_on_path$y_loc), color="yellow")+
  #geom_point(aes(x= max_locs_on_path$x_loc, y=-max_locs_on_path$y_loc), color="yellow")+
  #geom_point(aes(x = mean(as.numeric(E10_5@meta.data$x_loc[E10_5@meta.data$annotation == "Heart"])), y = -mean(as.numeric(E10_5@meta.data$y_loc[E10_5@meta.data$annotation == "Heart"]))), color = "#5D0E41", size = 2) +
  geom_point(data = windows, aes(x= sliding_x, y=sliding_y), color="#A0153E", stroke=6)+
  annotate("path",
           x=mean(as.numeric(E10_5@meta.data$x_loc[E10_5@meta.data$annotation == "Heart"]))+45*cos(seq(0,2*3.14,length.out=100)),
           y=-mean(as.numeric(E10_5@meta.data$y_loc[E10_5@meta.data$annotation == "Heart"]))+45*sin(seq(0,2*3.14,length.out=100)))+
  annotate("path",
           x=mean(as.numeric(E10_5@meta.data$x_loc[E10_5@meta.data$annotation == "Heart"]))+30*cos(seq(0,2*3.14,length.out=100)),
           y=-mean(as.numeric(E10_5@meta.data$y_loc[E10_5@meta.data$annotation == "Heart"]))+30*sin(seq(0,2*3.14,length.out=100))) 



bind_rows(sliding_ligand.df, .id = "source")->combined_df


# Calculate the mean and confidence intervals
summary_df <- combined_df %>%
  group_by(source) %>%
  summarise(
    mean_ligand = mean(ligand, na.rm = TRUE),
    lower_ci = mean_ligand - qt(0.975, df = n() - 1) * sd(ligand, na.rm = TRUE) / sqrt(n()),
    upper_ci = mean_ligand + qt(0.975, df = n() - 1) * sd(ligand, na.rm = TRUE) / sqrt(n())
  )

# Plot using ggplot2
ggplot(summary_df, aes(x = as.numeric(source), y = mean_ligand, group = 1)) +
  geom_line(color = "#E1AFD1") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), color="#AD88C6", alpha = 0.2) +
  labs(title = "Ptch1 - E10.5") +
  theme_classic()+
  #scale_color_manual(values = c("#E1AFD1", "#AD88C6"))+
  #geom_hline(yintercept =0)+
  theme(
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    #legend.position = "none",
    plot.title = element_text(hjust = 0.5, size=20),
    #axis.text = element_blank(),
    text = element_text(size=18))+labs(x="Location from SHF-OFT", y="Ptch1", tag="C")->plt.marker.C





plt.density.A + labs(x="Enrichment", y="Density") + ggtitle("Osr1/Ptch1") + plt.heart.vicinity+ labs(x="x", y="y") + plt.marker.C + labs(x="SHF-OFT")



write.csv(summary_df, "~/Downloads/E10_5_ptch1.csv")
















