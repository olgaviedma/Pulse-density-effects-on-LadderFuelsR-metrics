# Impact of LiDAR Pulse Density on Forest Fuels Metrics derived using LadderFuelsR

## Authors: Olga Viedma and JM Moreno

### National LiDAR surveys conducted in various countries often use low pulse densities, which can significantly impact the accuracy of derived fuel property estimates, potentially leading to underestimation or overestimation of critical metrics used for fire hazard assessments. This study evaluates the impact of LiDAR pulse density thinning on forest structure characterization using the LadderFuelsR metrics. High-resolution LiDAR data from some Mediterranean forests was systematically thinned to simulate varying pulse densities. Key forest metrics, including leaf area density (LAD), leaf area index (LAI), canopy base height, depth, height of fuel layers, and inter-layer distances, were derived to assess changes at different thinning levels.

## Getting Started

# SECTION 1. THINNING, CHM, SEGMENTATION, LAD PROFILES <br/>
### Download ZENODO database: DATA.zip. Available at: https://zenodo.org/record/14786024
```{r pressure, echo=FALSE}
# Load necessary library
library(httr)

# Step 1: Define URLs and Local Paths to download ZENODO database

zenodo_url <- "https://zenodo.org/record/14786024/files/DATA.zip?download=1"  # Zenodo download link
zip_file <- "DATA.zip"       # Name of the downloaded ZIP file
output_dir <- "DATA"         # Final directory for extracted files
temp_dir <- "temp_extracted" # Temporary folder for intermediate extraction

# Step 2: Download the ZIP File (if not already downloaded)

if (!file.exists(zip_file)) {
  message("Downloading DATA.zip from Zenodo...")
  GET(zenodo_url, write_disk(zip_file, overwrite = TRUE))
} else {
  message("DATA.zip already exists. Skipping download.")
}

# Step 3: Validate the Downloaded ZIP File

file_size <- file.info(zip_file)$size
if (is.na(file_size) || file_size < 1e6) {  # Ensure the file is >1 MB
  stop("Download failed: ZIP file is too small or empty!")
} else {
  message(paste("ZIP file size:", round(file_size / (1024 * 1024), 2), "MB"))
}

# Step 4: Extract the ZIP File to a Temporary Folder

message("Extracting the ZIP file to a temporary folder...")
unzip(zip_file, exdir = temp_dir, overwrite = TRUE)

# Step 5: Move Contents to Final Output Directory

if (file.exists(file.path(temp_dir, "DATA"))) {
  # If ZIP contains an extra "DATA" folder, move its contents
  message("Moving extracted contents to the DATA directory...")
  file.rename(from = file.path(temp_dir, "DATA"), to = output_dir)
} else {
  # Otherwise, directly move the extracted files
  message("Moving extracted files to the DATA directory...")
  file.rename(from = temp_dir, to = output_dir)
}

# Step 6: Clean Up Temporary Files

unlink(temp_dir, recursive = TRUE)  # Delete the temporary folder

# Step 7: Verify the Final Extracted Files

message("Veifying the final structure of extracted files...")
extracted_files <- list.files(output_dir, recursive = TRUE)
if (length(extracted_files) == 0) {
  stop("Unzipping failed: No files were extracted.")
} else {
  message("Extraction successful! Files extracted:")
  print(extracted_files)
}

# Step 8 (Optional): Print a Preview of Extracted Files

message("First few extracted files:")
print(head(extracted_files))

```
## 1.THINNING LIDARPOD FILES
### The LiDAR point clouds were thinned from the original density (median 648 pulses/m2 ± 376.4) to only 1 pulses/m2 by using these thinned levels (pulse/m2): 100, 50, 25, 10, 5, 4, 3, 2 and 1 in LasTools. We used the “random” option selecting the lowest (default) points.
## Full resolution vs 1 pulse/m²
![Full resolution vs 1 pulse/m²](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/1_THINNING_LIDARPOD.png)
## 2. CHM FOR LAS FILES
### The Canopy Height Model (CHM) was calculated using the pitfree method, scanning the canopy at different height intervals (0 to 40 meters) and using variable-sized windows at 0.5 m resolution 
## CHM (Full resolution vs 1 pulse/m²
![CHM (Full resolution vs 1 pulse/m²)](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/2_CHMs.png)
## 3-4. CROWNS DELIMITATION AND METRICS FOR UNTHINNED LAS FILES
### Individual tree crowns were derived using the watershed algorithm implemented in LidR package. The "minimum height"" was set at 4 meters. The "tolerance"" parameter was set at 0.25 and the parameter “ext” at 1.
## Crowns polygons at full resolution (watershed)
![Crowns polygons at full resolution (watershed)](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/3_CROWN_POLYGONS.png)
## 5. CROP ALL THINNED .LAS FILES WITH TREE POLYGONS DERIVED FROM BENCHMARK .LAS FILE (UNTHINNED)
### All returns from the normalized LiDAR heights were cropped using the crown’s polygons from the benchmark Lidar file (unthinned)
## Cropping all thinned LiDAR files with full-resoultion crowns polygons
![Cropping all thinned LiDAR files with full-resoultion crowns polygons](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/4_CROPPING_LAS.png)
## 6_7. LAI-LAD METRICS BY TREE
### The Leaf Area Density (LAD) profiles at 1-meter intervals (height bins) were calculated for each tree
## LAD profiles and LAI metrics
![LAD profiles and LAI metrics](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/5_LAD_LAI_METRICS.png)

# SECTION 2.LadderFuelsR and LAI metrics (STEPS 8-19)<br/>
## 18. JOINING LADDER FUEL PROPERTIES TO CROWN POLYGONS
### We applied the methodology developed in LadderFuelsR to obtain several relevant variables of the vertical structure of the trees (Viedma et al. 2024), following a sequential workflow (https://github.com/olgaviedma/LadderFuelsR) (accessed on 11th November 2024). 
## LadderFuelsR metrics associated to crown polygons
![LadderFuelsR metrics associated to crown polygons](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/18_ladderfuelsR_metrics_polygons.png)
## 19. ALL METRICS and DIFFERENCES: (NO THINNED - REST OF THINNED LEVELS)
### Internal tables for further processing

# SECTION 3. CLUSTERING ALL TREES (across thinning levels)<br/>
## 20. PCA CLUSTER AND FREQUENCY DISTRIBUTION (FIGS.2 & 4 and FIG.S3)
### The fuel properties of each tree at each thinned level (n=10) were classified using the Hierarchical Clustering on Principal Components (HCPC) from the FactoMineR package. The performance of different numbers of clusters was tested using: the silhouette scores and two global performance indices: the Calinski-Harabas and the Davies-Bouldin. To assess the susceptibility of the forest structures to change to other ones over thinned levels, we calculated the occupancy percentage of each cluster at each thinning level.
## Figure 2. Clusters on PCA axes
![Figure 2. Clusters on PCA axes](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/FIG2_plot_clusters_PCA_AXES.png)
## Table S2. Clusters Performance
![Table S2. Clusters performance](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/TABLE_S2_CLUSTERS_PERFORMANCE.PNG)
## Figure S3. Clusters distribution by sites
![Figure S3. Clusters distribution by sites](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/FIG_S3_plot_distribution_cluster_by_CLUST_NT_by_thin_levels.png)
## 21. STATISTICAL ANALYSIS "PCA CLUSTER": kruskal-wallis (FIG. S2  TABLE S3)
### The clusters statistical separability was checked using the non-parametric post hoc ANOVA test Kruskal–Wallis. 
## Figure S2. Statistical differences among clusters (Kruskal-Wallis test)
![Figure S2. Statistical differences among clusters (Kruskal-Wallis test)](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/FIG_S2_plot_CLUSTER_estadisticos_kruskal.png)
## Table S3. Statistics of main fuel properties by cluster
![Table S3. Statistics of main fuel properties by cluster](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/TABLE_S3_CLUSTERS_STATISTICS.PNG)
## 22. JOINING PCA CLUSTERS AND FITTED MODELS TO CROWN POLYGONS (input for FIG. S4)
## Figure S4. Assignation of clusters to crowns polygons across thinned levels
![Figure S4. Assignation of clusters to crowns polygons across thinned levels](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/FIG_S4_plot_clusters_on_crown_polygons.png)

# SECTION 4.CALCULATE MEDIAN VALUES
### Internal tables for further processing

# SECTION 5. REGRESSION MODELS USING MEDIAN VALUES DIRECTLY: RATES OF CHANGE (TABLE S4 & FIG. 7)
### The rates of change (RCH) of the median fuel properties (aggregated from all trees) across the thinned levels (n=10) were calculated using regression models (exponential, logarithmic, and linear). The resulting regression slopes were then expressed as percentages.
## Table S4. Rates of Change (RCH) of main fuel properties using median values directly
![Table S4. Rates of Change (RCH) of main fuel properties using median values directly](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/TABLE_S4_RCH_by_MEDIAN.PNG)
## Figure 7. Rates of Change (RCH) of main fuel properties using median values directly
![Figure 7. Rates of Change (RCH) of main fuel properties using median values directly](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/FIG7_PLOT_RATES_CHANGE_MEDIAN_VAR_p2575.png)

# SECTION 6. REGRESSION MODELS USING INDIVIDUAL TREES: RATES OF CHANGE RATES FOR EACH TREE (TABLE S5 & FIG S6)
### The rates of change (RCH) in fuel properties for each individual tree across the thinned levels were calculated using multiple regression models (exponential, logarithmic, and linear). The resulting regression slopes were then converted into percentage values.
## Table S5. Rates of Change (RCH) of main fuel properties from individual trees
![Table S5. Rates of Change (RCH) of main fuel properties from individual trees](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/TABLE_S5_RCH_by_TREES.PNG)
## Figure S6. Comparison of Rates of Change (RCH) of main fuel properties calculated from median values and from individual trees
![Figure S6. Comparison of Rates of Change (RCH) of main fuel properties calculated from median values and from individual trees](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/FIG_S6_PLOT_COMPARISON_RATES_CHANGE_TREES_MEDIAN_VAR_p2575.png)

# SECTION 7. RATES OF CHANGE BY CLUSTERS: TAKING THE "NO_THINNED" LEVEL CLASSIFICATION AS BENCHMARK
### To assess if certain tree typologies identified at the highest resolution underwent significant changes at other thinning levels, we calculated the percentage of trees that changed their cluster assignment at full resolution across the different thinning levels.
## Figure S3. Percentage of trees that changed their full resolution cluster assignment across thinning levels
![Figure S3. Percentage of trees that changed their full resolution cluster assignment](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/FIG_S3_plot_distribution_cluster_by_CLUST_NT_by_thin_levels.png)
## 32. MEDIAN VALUES CLUSTERS (CLUST_NT)
### Internal tables for further processing
## 33.RATES OF CHANGE BY CLUSTERS (USING THE BENCHMARK: CLUST_NT) (TABLE S6 & FIGS. S8-S12)
### To estimate the RCH at the cluster level, we tracked the trees belonging to the clusters identified at full resolution across the thinning levels. Following, regression models were applied to each tree individually by clusters to get the median RCH from the best model (with the highest adjusted R²).
## Table S6. Rates of Change (RCH) of main fuel properties by Cluster using median values
![Table S6. Rates of Change (RCH) of main fuel properties by Cluster using median values](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/TABLE_S6_RCH_by_CLUSTERS_NT_BY_MEDIANS.PNG)
## Figures S8-12. Rates of Change (RCH) of main fuel properties by Cluster using median value
![Figures S8-12. Rates of Change (RCH) of main fuel properties by Cluster using median value](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/FIG_S8_plot_RCH_quantiles_2575_CLUSTER_NT_F1.png)
## 34.COMPARISON OF THE RATES OF CHANGE FROM INDIVIDUAL TREES vs. MEDIAN VALUES BY CLUSTERS (TABLE S7 & FIG. S7)
## Table S7. Rates of Change (RCH) of main fuel properties by Cluster using all trees
![Table S7. Rates of Change (RCH) of main fuel properties by Cluster using all trees](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/TABLE_S7_RCH_by_CLUSTERS_NT_BY_TREES.PNG)
## Figures S7.Comparison of Rates of Change (RCH) "by Clusters" of main fuel properties calculated from median values and for each individual tree
![Figures S7. Comparison of Rates of Change (RCH) "by Clusters" of main fuel properties calculated from median values and for each individual tree](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/FIG_S7_PLOT_COMPARISON_RATES_CHANGE_CLUSTER_NT_TREES_MEDIAN_p2575.png)

## PLOT EFFECTIVE METRICS (LadderFuelsR) by CLUSTERS (inputs for FIG. 6)
## Figure 6. Some examples of changes in cluster assignment from full-resolution classification to the lowest thinning level (1 pulses/m2). 
![Figure 6. Some examples of changes in cluster assignment from full-resolution classification to the lowest thinning level (1 pulses/m2). ](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/Fig6_upper_panel.png)

# Acknowledgements

We gratefully acknowledge funding from project INFORICAM (PID2020-119402RB-I00), funded by the Spanish MCIN/AEI/ 10.13039/501100011033 and by the "European Union NextGenerationEU/PRTR". Carlos Silva was supported by the NASA's Carbon Monitoring System funding (CMS, grant 22-CMS22-0015).

# Reporting Issues

Please report any issue regarding this study to Dr. Olga Viedma ([olga.viedma\@uclm.es](mailto:olga.viedma@uclm.es){.email})

# Citing this paper

Viedma,O. & Moreno, JM: Impact of LiDAR Pulse Density on Forest Fuels Metrics derived using LadderFuelsR. Ecological Informatics

