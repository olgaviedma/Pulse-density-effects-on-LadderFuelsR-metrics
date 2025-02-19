
#**Impact of LiDAR Pulse Density on Forest Fuels Metrics derived using LadderFuelsR**

#**Authors: Olga Viedma and JM Moreno**

#This study evaluates the impact of LiDAR pulse density thinning on forest structure characterization using the LadderFuelsR framework. High-resolution LiDAR data were collected from Mediterranean forest sites with diverse structures and fire histories and systematically thinned to simulate varying pulse densities. Key forest metrics, including leaf area density (LAD), leaf area index (LAI), canopy base height, depth, height of fuel layers, and inter-layer distances, were derived to assess changes at different thinning levels.
#This third section deals with the classification of each tree at each thinning level (n= 10) using the fuel properties derived from LadderFuelsR. 
#The fuel properties at each thinned level were scaled by z-scores and classified using the Hierarchical Clustering on Principal Components (HCPC) from the FactoMineR package. 
#The performance of different numbers of clusters was checked using various tests: the silhouette scores and two global performance indices: the Calinski-Harabas (higher values, better the clustering) and the Davies-Bouldin (lower values, better the clustering). 
#The main statistics (median, min, max, and standard deviation) were extracted, and their statistical separability was tested using the non-parametric post hoc ANOVA test Kruskal–Wallis. 
#The occupancy percentage of each cluster at each thinning level was calculated. 

# Getting Started

## Installation

remotes::install_github("Jean-Romain/lidRplugins")

#The CRAN version:
install.packages("LadderFuelsR")

# The development version:
install.packages("remotes")
library(remotes)
install_github("https://github.com/olgaviedma/LadderFuelsR", dependencies = TRUE)

install.packages("rgdal", dependencies = TRUE)
library(remotes)
library(rgdal)
install_github("r-spatial/sf")


## Required libraries

library(devtools)
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)
library(stringi)
library(purrr)
library(rlang)
library(tidyverse)
library(magrittr)
library(rlist)

library(sp)
library(sf)
library(raster)
library(data.table)
library(rgdal)

library(lidR)
library(leafR)
library(lidRplugins)
library(LadderFuelsR)


########################################################
#20. PCA CLUSTER AND FREQUENCY DISTRIBUTION (FIGS.3 & 5 and TABLE S2 & INPUT FOR FIG.S4)
########################################################

library(tidyverse)
library(FactoMineR)
library(factoextra)
library(cluster)
library(clusterCrit)
library(fpc)
library(ggrepel)

library(officer)
library(flextable)


# INPUT DATA

cluster_dir<-dir.create("RESULTS/CLUSTER", recursive = TRUE )
cluster_dir<-"RESULTS/CLUSTER"

all_data_lidarpod1 <- read.table("WATER_LAI/ALL_METRICS_together_long_NO_THINNED_depurated.txt", sep = "\t", header = TRUE)


all_data_lidarpod1$treeID <- factor(all_data_lidarpod1$treeID,
                                    levels=unique(all_data_lidarpod1$treeID))
all_data_lidarpod1$thin_level <- factor(all_data_lidarpod1$thin_level,
                                        levels=unique(all_data_lidarpod1$thin_level))
all_data_lidarpod1$sites <- factor(all_data_lidarpod1$sites,
                                   levels=unique(all_data_lidarpod1$sites))


# Set xlim and ylim parameters outside the loop
max_x <- NULL
max_y <- NULL

# Get unique thin_levels excluding "PNOA"
unique_thin_levels <- unique(all_data_lidarpod1$thin_level)


all_data_lidarpod1b <- all_data_lidarpod1[, c("treeID","treeID1","thin_level" ,"sites", "Hcbh1","dptf1","Hcb1_H1","Hdptf1",  "mxld_Hc","bp_Hcbh","Hcbh_br","mxld_dp","bp_dptf","mxld_ff","bp_ffds", "mxld_Hdp", "bp_Hdpt","mx_hgh_", "nlayers", "lai_tot", "undrst_")]

names(all_data_lidarpod1b) <- c("treeID","treeID1","thin_level" ,"sites", "F1_cbh", "F1_depth", "F1_LAD","F1_Hdepth","MXLAD_cbh", "BP_cbh", "BR_cbh","MXLAD_depth", "BP_depth", "MXLAD_dist", "BP_dist","MXLAD_Hdepth","BP_Hdepth", "MAX_height", "Nº_layers", "total_LAI", "understory_LAI")


all_data_lidarpod1c <- all_data_lidarpod1b[, -c(1:4)]

# Display the updated dataframe
print(all_data_lidarpod1c)


## PCA CLUSTER


# Prepare data and scale it
all_data_lidarpod2 <- all_data_lidarpod1c %>%
  na.omit() %>%          # Remove missing values (NA)
  scale()                # Scale variables

# PCA
res.pca <- PCA(all_data_lidarpod2,  ncp = Inf)


#########3

# Create the PCA plot for variables only
pca_var_plot <- fviz_pca_var(
  res.pca,
  col.var = "black",   # Color of variable arrows and labels
  repel = TRUE,        # Avoid overlapping labels
  labelsize = 5        # Adjust the size of the variable labels here
) + 
  theme_minimal() +
  theme_bw() +
  theme(
    axis.text = element_text(size = 16),      # Increase size of axis tick labels
    axis.title = element_text(size = 18),     # Increase size of axis titles
    legend.text = element_text(size = 16),    # Increase size of legend text if used
    legend.title = element_text(size = 16),   # Increase size of legend title
    plot.title = element_text(size = 20),     # Increase size of plot title
    text = element_text(size = 16)            # Increase size of other text elements
  )

# Display the plot
pca_var_plot

# Define the file path
plot_file_path <-  file.path(cluster_dir, "PCA_plot.png")

# Save the plot
ggsave(
  filename = plot_file_path, 
  plot = pca_var_plot, 
  width = 8, 
  height = 6, 
  dpi = 300
)

#########################3

doc <- read_docx()

eigenvalues<- res.pca$eig

# Round the values to 2 decimal places
eigenvalues <- round(eigenvalues, 2)
eigenvalues <- as.data.frame(eigenvalues)
eigenvalues$component <- rownames(eigenvalues)
eigenvalues<- eigenvalues[,c(4,1:3)]

# Reset the row names
rownames(eigenvalues) <- NULL

ft <- flextable(eigenvalues)
autofit(ft)  # Automatically adjusts the column width to fit the content

# Add the flextable to the document
doc <- doc %>% 
  body_add_flextable(ft)


# Save the Word document
output_path <- file.path(cluster_dir,"PCAs_eigenvalues.docx")
print(doc, target = output_path)

##########3

doc1 <- read_docx()

variables_corr<- res.pca$var$cor

# Round the values to 2 decimal places
variables_corr <- round(variables_corr, 2)
variables_corr <- as.data.frame(variables_corr)
variables_corr$Variables <- rownames(variables_corr)
variables_corr<- variables_corr[,c(18,1:5)]
variables_corr <- variables_corr %>%
  arrange(desc(Dim.1))

# Reset the row names
rownames(variables_corr) <- NULL

ft1 <- flextable(variables_corr)
autofit(ft1)  # Automatically adjusts the column width to fit the content

# Add the flextable to the document
doc1 <- doc1 %>% 
  body_add_flextable(ft1)


# Save the Word document
output_path <- file.path(cluster_dir,"PCAs_var_corr.docx")
print(doc1, target = output_path)

#############################3

# Update max_x and max_y with maximum absolute values of principal components
max_x <- max(max_x, max(abs(res.pca$ind$coord[, 1])), na.rm = TRUE)
max_y <- max(max_y, max(abs(res.pca$ind$coord[, 2])), na.rm = TRUE)

# HIERARCHICAL ASCENDANT CLUSTERING
res.hcpc <- HCPC(res.pca, kk=Inf, 4, consol=TRUE)
clust_pca <- res.hcpc$data.clust
cluster_memberships <- res.hcpc$data.clust$clust

cluster_memberships_numeric <- as.numeric(as.character(cluster_memberships))


## INPUT FOR TABLE S2: V.TEST (VARIABLES IMPORTANCE)

# Load necessary packages
library(officer)
library(dplyr)

# Assuming res.hcpc$desc.var$quanti is your input
quanti_data <- res.hcpc$desc.var$quanti

# Initialize an empty list to store data for each cluster
all_data <- list()

# Loop through each cluster (e.g., `1`, `2`, etc.)
for (cluster in names(quanti_data)) {
  # Extract relevant data
  cluster_data <- quanti_data[[cluster]]
  
  # Convert matrix to data frame
  cluster_data_df <- as.data.frame(cluster_data)
  
  cluster_data_df <-cluster_data_df %>% mutate(across(where(is.numeric), ~ formatC(., format = "f", digits = 2))) 
  
  cluster_data_df$Variable = rownames(cluster_data_df)
  names(cluster_data_df$Variable) = rownames(cluster_data_df)
  
  # Create a data frame with the necessary columns
  cluster_table <- cluster_data_df %>%
    dplyr::select(Variable,  # Use the row names for the Variable column
                  v.test,
                  p.value) %>%
    mutate(Cluster = cluster)  # Add a column for cluster identification
  
  # Append the cluster table to the list
  all_data[[cluster]] <- cluster_table
}

# Combine all data frames into one
final_table <- bind_rows(all_data)

# Reset the row names
rownames(final_table) <- NULL

final_table <- final_table[,c(4,1:3)]

# Create a Word document
doc2 <- read_docx()

ft2 <- flextable(final_table)
autofit(ft2)  # Automatically adjusts the column width to fit the content

# Add the flextable to the document
doc2 <- doc2 %>% 
  body_add_flextable(ft2)

# Save the Word document
output_path <- file.path(cluster_dir,"TABLE_S2_PCAs_var_importance_vtest.docx")
print(doc2, target = output_path)


# INPUT FOR TABLE S2: CLUSTERS PERFORMANCE

library(clusterSim)

# Calculate silhouette scores
silhouette_scores <- silhouette(cluster_memberships_numeric, dist(all_data_lidarpod2))
# Extract silhouette scores for each cluster
cluster_silhouette_scores <- split(silhouette_scores[, "sil_width"], cluster_memberships_numeric)

cluster_silhouette_score2<-list()
for (cluster_id in unique(cluster_memberships_numeric)) {
  cluster_silhouette_score <- cluster_silhouette_scores[[as.character(cluster_id)]]
  cluster_silhouette_score1<-paste("Silhouette Score for Cluster", cluster_id, ":", mean(cluster_silhouette_score))
  cluster_silhouette_score2[[cluster_id]] <- cluster_silhouette_score1
}

# Split the character strings and extract numeric values
cluster_silhouette_score_numeric <- sapply(cluster_silhouette_score2, function(x) as.numeric(sub(".*:", "", x)))
cluster_silhouette_score3 <- as.data.frame(cluster_silhouette_score_numeric)
cluster_silhouette_score3$cluster <- paste0("Cluster_", seq_len(nrow(cluster_silhouette_score3)))

# Calinski-Harabasz Index: A higher value of CH indicates a better clustering, because it means that the data points are more spread out between clusters than they are within clusters.

calinski_harabasz_index <- cluster.stats(dist(all_data_lidarpod2), cluster_memberships_numeric)$ch
calinski_harabasz_index1 <- as.data.frame(calinski_harabasz_index)

db_index <- index.DB(all_data_lidarpod2, cluster_memberships_numeric)$DB
print(db_index)
db_index1 <- as.data.frame(db_index)

cluster_performance_all<- RbindAll(cluster_silhouette_score3,calinski_harabasz_index1,db_index1)

### Input for SM Table 2. 

write.table(cluster_performance_all,file.path(cluster_dir,"TABE_S2_performance_PCA4_all_metrics.txt"), row.names = F, sep = "\t")


# CLUSTERS STATISTICS

clust1_median <- aggregate(all_data_lidarpod1c, by=list(cluster=clust_pca$clust), median)
clust1_min <- aggregate(all_data_lidarpod1c, by=list(cluster=clust_pca$clust), min)
clust1_max <- aggregate(all_data_lidarpod1c, by=list(cluster=clust_pca$clust), max)
clust1_deviat <- aggregate(all_data_lidarpod1c, by=list(cluster=clust_pca$clust), sd)

clust_PCA_statist <- cbind(clust1_median, clust1_deviat, clust1_min, clust1_max)

clust_pca1 <- cbind(all_data_lidarpod1b, cluster = clust_pca$clust)

# Define the order of thin_level levels
new_order <- c("NO_THINNED",
               "THINNED_100P",
               "THINNED_50P",
               "THINNED_25P",
               "THINNED_10P",
               "THINNED_5P",
               "THINNED_4P",
               "THINNED_3P",
               "THINNED_2P",
               "THINNED_1P")

# Reorder factor levels
clust_pca1$thin_level <- factor(clust_pca1$thin_level, levels = new_order)

write.table(clust_pca1, file = file.path(cluster_dir, "ALL_METRICS_all_sites_CLUST_PCA4.txt"), row.names = F, sep = "\t")

write.table(clust_PCA_statist, file.path(cluster_dir,"statistics_PCA4_all_metrics.txt"), row.names = F, sep = "\t")



# FIG.3: PLOT CLUSTERS IN PCA AXES 


plot_file1 <- file.path(cluster_dir,"FIG3_plot_clusters_PCA_AXES.png")

cluster_colors <- c("1" = "orange", "2" = "brown","3" = "green","4" = "lightblue")  # Initialize the color vector

library(factoextra)  # Ensure the factoextra library is loaded for fviz_pca_biplot
library(ggplot2)     # Ensure ggplot2 is loaded for additional customizations
library(ggrepel)     # Load ggrepel for geom_text_repel

# Create a data frame for PCA variable coordinates
pca_var_coords <- as.data.frame(get_pca_var(res.pca)$coord)

plot1 <- fviz_pca_biplot(
  res.pca, 
  geom = "point", 
  pointsize = 1, 
  label = "var", 
  col.var = "black", 
  habillage = clust_pca1$cluster, 
  addEllipses = TRUE, 
  ellipse.level = 0.95, 
  repel = TRUE, 
  show.legend = FALSE,  
  labelsize = 4       # Increase label size
) +
  xlim(c(-11, 11)) + 
  ylim(c(-8, 8)) + 
  scale_color_manual(values = cluster_colors, breaks = names(cluster_colors)) + 
  theme_bw() +  
  theme(
    axis.text = element_text(size = 14),      # Increase tick label size
    axis.title = element_text(size = 14),     # Increase axis title size
    legend.text = element_text(size = 14),    # Increase legend text size
    legend.title = element_text(size = 14),   # Increase legend title size
    plot.title = element_text(size = 16),     # Increase plot title size
    text = element_text(size = 14)            # Increase overall text size
  ) +
  coord_equal()


# Print the plot
print(plot1)

# Save plot
ggsave(filename = plot_file1, plot = plot1, width = 8.68, height = 7.35)


##FIG.5 (CLUSTERS DISTRIBUTION BY THINNING LEVELS)

library(epiDisplay)
library(questionr)
library(datawizard)
library(sjPlot)
library(parameters)

clust_pca1 <- read.table(file.path(cluster_dir,"ALL_METRICS_all_sites_CLUST_PCA4.txt"), sep = "\t", header = TRUE)

# Reorder the thin_level levels to be in reverse order
clust_pca1$thin_level <- factor(clust_pca1$thin_level)

# Define the order of thin_level levels
new_order <- c("NO_THINNED",
               "THINNED_100P",
               "THINNED_50P",
               "THINNED_25P",
               "THINNED_10P",
               "THINNED_5P",
               "THINNED_4P",
               "THINNED_3P",
               "THINNED_2P",
               "THINNED_1P")

clust_pca1$thin_level <- factor(clust_pca1$thin_level, levels = new_order)
clust_pca1 <- clust_pca1[order(clust_pca1$thin_level), ]

new_names <- c(">300",
               "100P",
               "50P",
               "25P",
               "10P",
               "5P",
               "4P",
               "3P",
               "2P",
               "1P")

# Rename levels
levels(clust_pca1$thin_level) <- new_names
clust_pca1 <- clust_pca1[order(clust_pca1$thin_level), ]



## FIG.4: PLOT OF CLUSTERS DISTRIBUTION (ALL SITES TOGETHER)  


plot_file2 <- file.path(cluster_dir,"FIG5_CLUSTER_distribution_ALL_sITES_TOGETHER.png")

# Create a ggplot object using ggplot2 directly
p1 <- ggplot(clust_pca1, aes(x = thin_level, fill = factor(cluster))) + 
  geom_bar(position = "fill") +  # Position "fill" creates stacked bar plot with proportional fill
  #coord_flip() +  # Flip coordinates for horizontal bars
  scale_fill_manual(values = c("orange", "brown", "green", "lightblue")) +  # Define custom colors for clusters
  labs(x = "Thinning levels", y = "Proportion", fill = "Cluster") +  # Set axis labels and legend title
  theme_minimal() +  # Use a minimal theme
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 16,color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 16,color = "black"),
        legend.title = element_text(size = 16,color = "black"),
        legend.text = element_text(size = 16,color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5,color = "black"),
        strip.text = element_text(size = 16,color = "black"))  # Increase the size of facet labels

# Display the plot
print(p1)

# Save the plot to file
ggsave(filename = plot_file2, plot = p1,width = 10, height = 8, dpi = 300)




## FREQUENCY STATISTICS OF CLUSTERS DISTRIBUTION (ALL SITES TOGETHER)


freq_table_clusters <- table(clust_pca1$thin_level,clust_pca1$cluster)
sum_freq <- apply(freq_table_clusters, 2, sum)
percentage_table <- data.frame(prop.table(freq_table_clusters, margin = 2) * 100)

names(percentage_table) <- c("thin_level", "cluster", "Freq_%")
write.table( percentage_table, file.path(cluster_dir, paste0("freq_distribution_cluster_PCA4_all_sites_together.txt")), sep = "\t", row.names = F)



## FIG. S4: PLOT OF CLUSTERS DISTRIBUTION FACET BY SITES 


clust_pca1$thin_level <- factor(clust_pca1$thin_level, levels = rev(levels(clust_pca1$thin_level)))

plot_file1 <- file.path(cluster_dir,"FIG_S4_freq_distribution_cluster_PCA4_separating_sites.png")

# Create a ggplot object using ggplot2 directly
p1 <- ggplot(clust_pca1, aes(x = thin_level, fill = factor(cluster))) + 
  geom_bar(position = "fill") +  # Position "fill" creates stacked bar plot with proportional fill
  coord_flip() +  # Flip coordinates for horizontal bars
  facet_wrap(~ sites, ncol = 2) +  # Create a panel for each vuelo
  scale_fill_manual(values = c("orange", "brown", "green", "lightblue")) +  # Define custom colors for clusters
  labs(x = "thin_level", y = "Proportion", fill = "Cluster") +  # Set axis labels and legend title
  theme_minimal() +  # Use a minimal theme
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 16,color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 16,color = "black"),
        legend.title = element_text(size = 16,color = "black"),
        legend.text = element_text(size = 16,color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5,color = "black"),
        strip.text = element_text(size = 16,color = "black"))  # Increase the size of facet labels

# Display the plot
print(p1)

# Save the plot to file
ggsave(filename = plot_file1, plot = p1,width = 10, height = 8, dpi = 300)



## FREQUENCY STATISTICS OF CLUSTERS DISTRIBUTION  BY SITES


# Create a frequency table grouped by thin_level, cluster, and sites
freq_table_clusters_sites <- table(clust_pca1$thin_level, clust_pca1$cluster, clust_pca1$sites)

# Calculate the sum of frequencies for each cluster within each site
sum_freq_sites <- apply(freq_table_clusters_sites, c(2, 3), sum)

# Calculate the percentage table for each site
percentage_table_sites <- as.data.frame(prop.table(freq_table_clusters_sites, margin = c(2, 3)) * 100)

# Rename columns for better clarity
colnames(percentage_table_sites) <- c("thin_level", "cluster", "sites", "percentage")

# Display the resulting table
head(percentage_table_sites)

write.table(percentage_table_sites, file.path(cluster_dir, paste0("freq_distribution_cluster_PCA4_separtaing_by_sites.txt")), sep = "\t", row.names = F)



##  OTHER RESULTS DERIVED FROM CLUSTER PCA  


options(scipen = 10)  # Adjust as needed
options(digits = 3)   # Adjust as needed

var_impor<-res.hcpc$desc.var
str(var_impor)
var_impor_list<-var_impor$quanti

clusters <-NULL
result<-NULL
for(i in seq_along(var_impor_list)){
  clusters<-do.call(rbind, var_impor_list[i])
  result<-rbind(result,clusters) 
}
result <- data.frame(result)
result <- result %>%
  dplyr::mutate_if(is.numeric, round,5)


write.table(result, file.path(cluster_dir,"var_importance_CLUST_PCA4.txt"), row.names=T, sep="\t")

####################3
pca_impor<-res.hcpc$desc.axes
str(pca_impor)
pca_impor_list<-pca_impor$quanti

clusters1 <-NULL
result1<-NULL
for(i in seq_along(pca_impor_list)){
  clusters1<-do.call(rbind, pca_impor_list[i])
  result1<-rbind(result1,clusters1) 
}

result1 <- data.frame(result1)
result1 <- result1 %>%
  dplyr::mutate_if(is.numeric, round,3)

write.table(result1, file.path(cluster_dir,"pca_importance_CLUST_PCA4.txt"), row.names=T, sep="\t")

####################3
pca_cor<-res.hcpc$call$t$res
corr_var<-pca_cor$var$cor

corr_var <- data.frame(corr_var)
corr_var <- corr_var %>%
  dplyr::mutate_if(is.numeric, round,3)

write.table(corr_var, file.path(cluster_dir,"var_corr_pca_CLUST_PCA4.txt"), row.names=T, sep="\t")


########################################################
#21.STATISTICAL ANALYSIS "PCA CLUSTER": kruskal-wallis (FIG. S2 & TABLE S3)
########################################################

library(car)
library(psych)
library(multcompView)
library(lsmeans)
library(FSA)
library(ggplot2)
library(phia)
library(multcomp)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(PMCMRplus)

# Turn off scientific notation globally
options(scipen = 999)

cluster_dir <- "RESULTS/CLUSTER"

all_data_lidarpod <- read.table(file.path(cluster_dir,"ALL_METRICS_all_sites_CLUST_PCA4.txt"), sep = "\t", header = TRUE)

# Reorder the thin_level levels to be in reverse order
all_data_lidarpod$thin_level <- factor(all_data_lidarpod$thin_level)

# Define the order of thin_level levels
new_order <- c("NO_THINNED",
               "THINNED_100P",
               "THINNED_50P",
               "THINNED_25P",
               "THINNED_10P",
               "THINNED_5P",
               "THINNED_4P",
               "THINNED_3P",
               "THINNED_2P",
               "THINNED_1P")

all_data_lidarpod$thin_level <- factor(all_data_lidarpod$thin_level, levels = new_order)
all_data_lidarpod <- all_data_lidarpod[order(all_data_lidarpod$thin_level), ]

new_names <- c(">300",
               "100P",
               "50P",
               "25P",
               "10P",
               "5P",
               "4P",
               "3P",
               "2P",
               "1P")

# Rename levels
levels(all_data_lidarpod$thin_level) <- new_names
all_data_lidarpod <- all_data_lidarpod[order(all_data_lidarpod$thin_level), ]


# Check the levels after reordering
levels(all_data_lidarpod$thin_level)
table(all_data_lidarpod$thin_level)

levels(all_data_lidarpod$sites)
table(all_data_lidarpod$sites)



# STATISTICS by CLUSTER


all_data_lidarpod1 <-all_data_lidarpod

all_data_lidarpod1$treeID<-factor(all_data_lidarpod1$treeID)
all_data_lidarpod1$treeID1<-factor(all_data_lidarpod1$treeID1)
all_data_lidarpod1$thin_level<-factor(all_data_lidarpod1$thin_level)
all_data_lidarpod1$cluster<-factor(all_data_lidarpod1$cluster)
all_data_lidarpod1$sites<-factor(all_data_lidarpod1$sites)


# Filter out non-numeric columns before summarizing
numeric_data <- all_data_lidarpod %>%
  dplyr::select( cluster, where(is.numeric)) %>%
  dplyr::select(-treeID1) # Exclude the treeID1 column


# Compute mean, standard deviation, and confidence intervals
median_data <- numeric_data %>%
  group_by(cluster) %>%
  summarize(
    across(
      everything(),
      list(
        median = ~median(.x, na.rm = TRUE),
        mean = ~mean(.x, na.rm = TRUE),
        min = ~min(.x, na.rm = TRUE),
        max = ~max(.x, na.rm = TRUE),
        sd = ~sd(.x, na.rm = TRUE),
        lower_ci = ~quantile(.x, 0.05, na.rm = TRUE),
        upper_ci = ~quantile(.x, 0.95, na.rm = TRUE)
      )
    ),
    .groups = 'drop'
  )



# Selecting only the median and related columns
median_data1 <- median_data %>%
  dplyr::select(ends_with(c("cluster", "_median", "_sd", "_min", "_max")))

# Check the column names to ensure the selection was correct
print(colnames(median_data1))

# Pivoting longer while keeping the "cluster" column separate
median_data_long <- median_data1 %>%
  pivot_longer(
    cols = -cluster,  # Exclude the "cluster" column from pivoting
    names_to = c("Variable", "ValueType"),  # Split variable names and types
    names_pattern = "(.*)_(median|sd|min|max)",  # Corrected the pattern
    values_to = "Value"  # Name for values column
  )

# Display the result
print(median_data_long)

median_data_long<-na.omit(median_data_long)

median_data_cluster_wide <- median_data_long %>%
  pivot_wider(
    names_from = ValueType,  # Move 'median' and 'sd' to new columns
    values_from = Value      # Values from the 'Value' column
  )

# Check the result
print(median_data_cluster_wide)

write.table (median_data_cluster_wide, file.path(cluster_dir, "MEDIAN_CLUSTERS_VARIABLES_QUANTILES0595_all_sites.txt"), sep = "\t", row.names =  F)



#KRUSKAL-WALLIS & ETA SQUARED by CLUSTER


library(dplyr)
library(tidyr)
library(multcompView)
library(ggplot2)

#Define the desired order of responses
responses <-c("F1_cbh", "F1_depth", "F1_LAD","F1_Hdepth","MXLAD_cbh", "BP_cbh", "BR_cbh",
              "MXLAD_depth", "BP_depth", "MXLAD_dist", "BP_dist",
              "MXLAD_Hdepth","BP_Hdepth", "MAX_height", "Nº_layers", "total_LAI", "understory_LAI")

results <- list()
eta_results <- list()
letters_results <- list()

for (k in seq_along(responses)) {
  cat("Processing response:", responses[k], "\n")
  
  # Check if all values are zero, skip this response if true
  if (all(all_data_lidarpod1[[responses[k]]] == 0)) {
    cat("Skipping", responses[k], "as all values are zero.\n")
    next
  }
  
  # Check if the response is a valid column
  if (responses[k] %in% colnames(all_data_lidarpod1)) {
    
    # Perform Kruskal-Wallis test
    ans <- kwAllPairsDunnTest(as.formula(paste(responses[k], "~ cluster")), data = all_data_lidarpod1)
    p_values <- ans$p.value
    
    # Filter out NA p-values and retain the corresponding pairs
    valid_comparisons <- which(!is.na(as.vector(p_values)))
    named_p_values <- as.vector(p_values)[valid_comparisons]
    
    # Generate names for the valid comparisons (e.g., "1-2", "1-3", ...)
    comparison_names <- apply(expand.grid(rownames(p_values), colnames(p_values)), 1, function(x) paste(sort(x), collapse = "-"))
    comparison_names <- comparison_names[valid_comparisons]
    
    # Assign names to the named_p_values vector
    names(named_p_values) <- comparison_names
    
    # Generate significance letters using multcompView
    sig_letters <- multcompLetters(p.adjust(named_p_values, method = "bonferroni"))
    
    
    # Convert the letters into a data frame
    letter_df <- data.frame(cluster = names(sig_letters$Letters), 
                            letter = sig_letters$Letters, 
                            variable = responses[k])
    letters_results[[k]] <- letter_df
    
    # Transform p-values to a long format for further processing or inspection
    p_values_df <- as.data.frame(p_values)
    p_values_df$cluster <- rownames(p_values_df)
    p_values_long <- tidyr::gather(p_values_df, comparison, p_value, -cluster, na.rm = TRUE)
    p_values_long$variable <- responses[k]
    p_values_long <- p_values_long %>% 
      mutate(across(where(is.numeric), ~ formatC(., format = "f", digits = 2)))
    
    results[[k]] <- p_values_long
    
    # Perform ANOVA and calculate eta squared
    model <- aov(as.formula(paste(responses[k], "~ cluster")), data = all_data_lidarpod1)
    eta1 <- as.data.frame(etaSquared(model, type = 2))
    eta1$variable <- responses[k]
    eta1 <- eta1 %>% 
      mutate(across(where(is.numeric), ~ formatC(., format = "f", digits = 2)))
    
    eta_results[[k]] <- eta1
    
  } else {
    cat("Column", responses[k], "not found in cluster_file1.\n")
  }
}

# Combine results for the current cluster
if (length(results) > 0) {
  results_kruskal <- do.call(rbind, results)
}
if (length(eta_results) > 0) {
  results_eta <- do.call(rbind, eta_results)
}
if (length(letters_results) > 0) {
  results_letters <- do.call(rbind, letters_results)
}

# Write final results to files
write.table(results_kruskal, file.path(cluster_dir, "kruskal_wallis_all_variables_by_CLUSTER_PCA4_allsites.txt"), sep = "\t", row.names = FALSE)
write.table(results_eta, file.path(cluster_dir,"eta_squared_all_variables_by_CLUSTER_PCA4_allsites.txt"), sep = "\t", row.names = FALSE)
write.table(results_letters,file.path(cluster_dir, "letters_kruskal_wallis_all_variables_by_CLUSTER_PCA4_allsites.txt"), sep = "\t", row.names = FALSE)


# TABLE S3 .DOCX STATISTICS CLUSTER: KRUSKAL-WALLIS 

library(officer)
library(flextable)

letter_kruskal <- read.table(file.path(cluster_dir, "letters_kruskal_wallis_all_variables_by_CLUSTER_PCA4_allsites.txt"), sep = "\t", header = TRUE)

median_data_cluster_wide <- read.table(file.path(cluster_dir, "MEDIAN_CLUSTERS_VARIABLES_QUANTILES0595_all_sites.txt"), sep = "\t", header = TRUE)

letter_kruskal$Variable <- as.character(letter_kruskal$variable)
median_data_cluster_wide$Variable <- as.character(median_data_cluster_wide$Variable)

letter_kruskal$cluster <- as.character(letter_kruskal$cluster)
median_data_cluster_wide$cluster <- as.character(median_data_cluster_wide$cluster)

# Ensure the Variable column is a factor with the desired levels
variables <- c( "F1_cbh","F1_depth","F1_LAD", "F1_Hdepth",  "MXLAD_cbh", "BP_cbh", "BR_cbh", 
                "MXLAD_depth", "BP_depth", 
                "MXLAD_dist", "BP_dist", 
                "MXLAD_Hdepth", "BP_Hdepth", "MAX_height", 
                "Nº_layers", "total_LAI", "understory_LAI"
)


median_data_cluster_wide$Variable <- factor(median_data_cluster_wide$Variable, levels = variables)
letter_kruskal$Variable <- factor(letter_kruskal$Variable, levels = variables)


# Order the dataframe by the factor levels of Variable
median_data_cluster_wide <- median_data_cluster_wide[order(median_data_cluster_wide$Variable), ]


# Perform the join
data_long2a <- dplyr::left_join(median_data_cluster_wide, letter_kruskal, by = c("cluster", "Variable"))
data_long2a <-data_long2a [,c(1:4,7,5:6)]

cluster_statistics_tab <- data_long2a %>% 
  mutate(across(where(is.numeric), ~ formatC(., format = "f", digits = 2)))

cluster_statistics_tab1 <- cluster_statistics_tab %>%
  dplyr::mutate(median = paste0(median, " ± ", sd),
                min_max = paste0("(", min, ",", max, ")")) %>%  # closed parentheses correctly
  arrange(Variable) %>%
  dplyr::select(cluster,Variable, median, letter, min_max)

# Reset the row names
rownames(cluster_statistics_tab1) <- NULL

# Create a Word document
doc3 <- read_docx()

ft3 <- flextable(cluster_statistics_tab1)
autofit(ft3)  # Automatically adjusts the column width to fit the content

# Add the flextable to the document
doc3 <- doc3 %>% 
  body_add_flextable(ft3)

# Save the Word document
output_path <- file.path(cluster_dir,"TABLE_S3_CLUSTER_estadisticos_kruskal.docx")
print(doc3, target = output_path)


## FIG. S2: PLOT CLUSTERS STATISTICS + KRUSKAL WALLIS LETTERS

library(ggplot2)

plot_file0a <- file.path(cluster_dir,"FIG_S2_plot_CLUSTER_estadisticos_kruskal.png")

# Initialize the color vector
cluster_colors <- c("1" = "orange", "2" = "brown", "3" = "green", "4" = "lightblue")

library(ggplot2)
library(ggrepel)  # Make sure to load ggrepel

plot0a <- ggplot(data_long2a, aes(x = as.factor(cluster), y = median, fill = as.factor(cluster))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = median - sd, ymax = median + sd), width = 0.2, position = position_dodge(0.9)) +
  geom_text(aes(y = median + sd + 0.1, label = letter),  # Increase offset to 0.1
            position = position_dodge(width = 0.9), 
            vjust = 0, 
            size = 4, 
            color = "black") +
  facet_wrap(~ Variable, scales = "free_y") +
  labs(title = "", 
       x = "Cluster", 
       y = "Mean Value ± SD", 
       fill = "Cluster") +
  scale_fill_manual(values = cluster_colors, breaks = names(cluster_colors)) +
  scale_y_continuous(limits = function(y) c(min(y) - 0.1 * diff(range(y)), max(y) + 0.2 * diff(range(y)))) +  # Extend y-axis limits
  theme_minimal() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5),
        strip.text = element_text(size = 14, face = "bold")) 

print(plot0a)

ggsave(plot0a, file = plot_file0a, width = 12, height = 8, dpi = 300)


############################################################
#22. JOINING PCA CLUSTERS AND FITTED MODELS TO CROWN POLYGONS (INPUTS FOR FIG.S4)
############################################################
library(dplyr)
library(sf)
library(stringr)

# Directory paths
cluster_dir <- "RESULTS/CLUSTER"
metrics_file_selec2 <- read.table(
  file.path(cluster_dir, "ALL_METRICS_all_sites_CLUST_PCA4.txt"),
  sep = "\t", header = TRUE
)

# Get the unique site names from the metrics data
sites_list <- unique(metrics_file_selec2$sites)

crown_dir <- "WATER_LAI"
# Get all crown subfolders, excluding "DATA_NO_THINNED"
metrics_crowns_subfolders <- list.dirs(crown_dir, full.names = TRUE, recursive = FALSE)
metrics_crowns_subfolders <- metrics_crowns_subfolders[!grepl("DATA_NO_THINNED", metrics_crowns_subfolders)]
metrics_crowns_subfolders1 <- list.dirs(metrics_crowns_subfolders, full.names = TRUE, recursive = FALSE)

# Loop over each site
for (site in sites_list) {
  cat("Processing site:", site, "\n")
  
  # Convert `sites` to a factor to ensure consistency
  metrics_file_selec2$sites <- factor(metrics_file_selec2$sites)
  
  # Filter metrics data for the current site (DataFrame)
  metrics_file_site <- subset(metrics_file_selec2, sites == site)
  
  # Ensure consistent column types
  metrics_file_site$treeID <- as.character(metrics_file_site$treeID)
  metrics_file_site$treeID1 <- as.character(metrics_file_site$treeID1)
  metrics_file_site$thin_level <- as.character(metrics_file_site$thin_level)
  
  
  # Loop over each subfolder to find shapefiles related to the current site
  for (subfolder in metrics_crowns_subfolders1) {
    # List all shapefiles in the current subfolder
    crown_files <- list.files(
      subfolder,
      pattern = glob2rx("ladderfuels_metrics_*_diff.shp"),
      full.names = TRUE,
      ignore.case = TRUE
    )
    
    if (length(crown_files) == 0) next  # Skip if no shapefiles are found
    
    for (crown_file in crown_files) {
      cat("Processing file:", crown_file, "\n")
      
      # Read the crown metric shapefile
      crown_metrics <- st_read(crown_file, quiet = TRUE)
      
      # Debug: Check shapefile structure
      cat("Shapefile structure:\n")
      print(summary(crown_metrics))
      
      # Add sites field if missing or incorrect
      if (!"sites" %in% colnames(crown_metrics)) {
        crown_metrics$sites <- gsub(".*_(SITE[0-9]+).*", "\\1", crown_file)
        cat("Inferred site from filename:", crown_metrics$sites[1], "\n")
      }
      
      # Ensure sites is a factor and filter for the current site
      crown_metrics$sites <- factor(crown_metrics$sites)
      crown_metrics <- subset(crown_metrics, sites == site)
      
      if (nrow(crown_metrics) == 0) {
        cat("No data for site:", site, "in file:", crown_file, "\n")
        next
      }
      
      # Validate geometry
      if (!all(st_is_valid(crown_metrics))) {
        crown_metrics <- st_make_valid(crown_metrics)
      }
      
      # Perform the join
      joined_data <- merge(
        crown_metrics, metrics_file_site,
        by = "treeID",
        all.x = TRUE
      )
      
      if (nrow(joined_data) == 0) {
        cat("No joined data for site:", site, "in file:", crown_file, "\n")
        next
      }
      
      # Drop redundant columns
      cols_to_drop <- c("treeID1.y", "treeIDf_y", "sites.y", "thin_level.y")
      joined_data <- joined_data[, !names(joined_data) %in% cols_to_drop]
      names(joined_data) <- sub("\\.x$", "", names(joined_data))
      
      # Convert to sf object
      joined_data_sf <- st_as_sf(joined_data, crs = st_crs(crown_metrics))
      
      # Validate geometry
      if (!all(st_is_valid(joined_data_sf))) {
        joined_data_sf <- st_make_valid(joined_data_sf)
      }
      
      # Function to abbreviate field names
      abbreviate_field_names <- function(field_names, max_length = 10) {
        make.unique(substr(field_names, 1, max_length))
      }
      
      # Apply abbreviation to field names
      names(joined_data_sf) <- abbreviate_field_names(names(joined_data_sf))
      
      # Save the output to a new shapefile
      output_file <- gsub(".shp$", paste0( "_CLUSTER4.shp"), crown_file)
      st_write(joined_data_sf, output_file, quiet = TRUE, append=F)
      cat("Saved joined file to:", output_file, "\n")
    }
  }
}
