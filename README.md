**Impact of LiDAR Pulse Density on Forest Fuels Metrics derived using LadderFuelsR**

Authors: Olga Viedma and JM Moreno

This study evaluates the impact of LiDAR pulse density thinning on forest structure characterization using the LadderFuelsR metrics.High-resolution LiDAR data collected from Mediterranean forest sites was systematically thinned to simulate varying pulse densities. 
Key forest metrics, including leaf area density (LAD), leaf area index (LAI), canopy base height, depth, height of fuel layers, and inter-layer distances, were derived to assess changes at different thinning levels.

# Getting Started

## SECTION 1. THINNING, CHM, SEGMENTATION, LAD PROFILES <br/>

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
### The LiDAR point clouds were thinned from the original density (median 648 pulses/m2 ± 376.4) to only 1 pulses/m2 by using these thinned levels (pulse/m2): 100, 50, 25, 10, 5, 4, 3, 2 and 1 in LasTools (REF). We used the “random” option selecting the lowest (default) points.
## Full resolution vs 1 pulse/m²
![Full resolution vs 1 pulse/m²](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/1_THINNING_LIDARPOD.png)

## 2. CHM FOR LAS FILES
### The Canopy Height Model (CHM) was calculated using the pitfree method, scanning the canopy at different height intervals (0 to 30 meters) and using variable-sized windows at 0.5 m resolution 
## CHM (Full resolution vs 1 pulse/m²
![CHM (Full resolution vs 1 pulse/m²)](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/2_CHMs.png)

## 3-4. CROWNS WATERSHED AND CROWNS METRICS FOR UNTHINNED LAS FILES
### Individual tree crowns was derived using the watershed algorithm implemented in LidR package. The minimum height was set at 4 meters. The tolerance parameter was set at 0.25 and the parameter “ext” at 1.
## [Crowns polygons at full resolution (watershed)
![Crowns polygons at full resolution (watershed)](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/3_CROWN_POLYGONS.png)

## 5. CROP ALL THINNED LAS FILES WITH TREE POLYGONS derived from BENCHMARK LAS FILE (UNTHINNED)
### All returns from the normalized LiDAR heights were cropped using the crown’s polygons from the benchmark Lidar file (unthinned)
## Cropping all thinned LiDAR files with full-resoultion crowns polygons
![Cropping all thinned LiDAR files with full-resoultion crowns polygons](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/4_CROPPING_LAS.png)

#6_7. LAI-LAD METRICS BY TREE
### The Leaf Area Density (LAD) profiles at 1-meter intervals (height bins) were calculated for each tree
## LAD profiles and LAI metrics
![LAD profiles and LAI metrics](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/5_LAD_LAI_METRICS.png)

#SECTION 2.LadderFuelsR METRICS AND LAI METRICS (STEPS 8-19)<br/>
#18. JOINING LADDER FUEL PROPERTIES WITH CROWN POLYGONS
```{r Joining crown polygons and ladder fuels metrics, echo=TRUE, message=FALSE, warning=FALSE}
```
# LadderFuelsR metrics associated to crown polygons
![LadderFuelsR metrics associated to crown polygons](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/18_ladderfuelsR_metrics_polygons.png)

#19.ALL METRICS and DIFFERENCES II: (NO THINNED - REST)
```{r Joining crown polygons and ladder fuels metrics, echo=TRUE, message=FALSE, warning=FALSE}
## Internal tables for further processing
```

#SECTION 3. CLUSTER OF ALL TREES (across thinning levels)<br/>
#20. PCA CLUSTER AND FREQUENCY DISTRIBUTION (FIGS.2 & 4 and FIG.S3)

```{r pressure, echo=FALSE}
```
# Table S2. Clusters Performance
![Table S2. Clusters performance](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/TABLE_S2_CLUSTERS_PERFORMANCE.PNG)
# Figure 2. Clusters on PCA axes
![Figure 2. Clusters on PCA axes](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/FIG2_plot_clusters_PCA_AXES.png)
# Figure S3. Clusters distribution by sites
![Figure S3. Clusters distribution by sites](https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/FIG_S3_plot_distribution_cluster_by_CLUST_NT_by_thin_levels.png)


#21. STATISTICAL ANALYSIS "PCA CLUSTER": kruskal-wallis (TABLE S3 & FIG. S2)

```{r pressure, echo=FALSE}

##############################################
# TABLE S3 .DOCX STATISTICS CLUSTER: KRUSKAL-WALLIS 
##############################################   

######################################
## FIG. S2: PLOT CLUSTERS STATISTICS + KRUSKAL WALLIS LETTERS
########################################
 
```

#22. JOINING PCA CLUSTERS AND FITTED MODELS TO CROWN POLYGONS (inputs for FIG. S4)

```{r Joining crown polygons and ladder fuels metrics, echo=TRUE, message=FALSE, warning=FALSE}

######################################
## FIG. S4: PLOT crowns polygons with CLUSTER assignation
########################################

```

#SECTION 4.CALCULATE MEDIAN VALUES

```{r pressure, echo=FALSE}

#23.CALCULATE MEDIAN I VALUES OF EACH VARIABLE BY TREES
#24. CALCULATE MEDIAN II VALUES OF EACH VARIABLE DIRECTLY

## Internal tables for further processing

```

#SECTION 5. REGRESSION MODELS USING MEDIAN VALUES DIRECTLY: RATES OF CHANGE (TABLE S4 & FIG. 7)

```{r pressure, echo=FALSE}

########################################
# TABLE S4: .DOCX BEST MODELS AND RATES OF CHANGE (SLOPE) USING MEDIAN VALUES
########################################

########################################
## FIG 7:PLOT RATES OF CHANGE BASED ON FITTED MODELS OVER MEDIAN VALUES
#############################################

```

#SECTION 6. REGRESSION MODELS USING INDIVIDUAL TREES: RATES OF CHANGE RATES FOR EACH TREE (median) (TABLE_S5 & FIG S6)

```{r, echo=TRUE, message=FALSE, warning=FALSE}



###########################################
#  TABLE S5. .DOCX BEST MODELS AND RATES OF CHANGE (SLOPE) USING ALL TREES
###########################################
  
###########################################
# FIG. S6: PLOTS OF COMPARISON BETWEEN RATES OF CHANGE (BAR PLOTS) FROM MODELS OF EACH ALL TREES AND FROM FITTED MODELS OVER MEDIANS DIRECTLY
 ################################


```

#SECTION 7. RATES OF CHANGE BY CLUSTERS: TAKING THE "NO_THINNED" LEVEL CLASSIFICATION AS BENCHMARK

```{r pressure, echo=FALSE}

#####################################
## FIG. S3.CLASS REASSIGNATION OF THE BENCHMARK CLUSTER (CLUSTER_NO THINNED) OVER THINNED LEVELS
#####################################

```

#32. MEDIAN VALUES CLUSTERS (CLUST_NT)

```{r pressure, echo=FALSE}


```

#33.RATES OF CHANGE BY CLUSTERS (USING THE BENCHMARK: CLUST_NT) (TABLE_S6 & FIGS. S8-S12)

```{r pressure, echo=FALSE}


###########################################
# TABLE_S6: .DOCX BEST MODELS AND RATES OF CHANGE (SLOPE) BY CLUSTER (USING MEDIAN VALUES) 
###########################################

######################################################################
# FIGS. S8-S12..PLOT FITTED MODELS AND RATES OF CHANGES BY CLUSTER (USING MEDIAN VALUES) 
######################################################################

 
```

#33.COMPARISON RATES OF CHANGE OF INDIVIDUAL TREES vs. MEDIAN VALUES (CLUST_NT) (TABLE S7 & FIG. S7)

```{r pressure, echo=FALSE}


###########################################
#  TABLE S7. CONVERT BEST MODELS FROM ALL TREES (BY thin_level + CLUSTER ) AND RATES OF CHANGE (SLOPE) INTO A .DOCX TABLE
###########################################

######################################
## FIG. S7: PLOT COMPARISIONS OF RATES OF CHANGE BY CLUSTER_NT + KRUSKAL WALLIS  (all trees vs median variables)
########################################


```

#00. PLOTS EFFECTIVE METRICS (LadderFuelsR) by CLUSTERS (inputs for FIG. 6)

```{r Plots of fuel layers with LAD percentage greater than 25 and the canopy base height (CBH) based on the maximum LAD percentage, echo=TRUE, message=FALSE, warning=FALSE}

```
