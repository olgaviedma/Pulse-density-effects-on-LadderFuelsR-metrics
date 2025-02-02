

**Impact of LiDAR Pulse Density on Forest Fuels Metrics derived using LadderFuelsR**

Authors: Olga Viedma and JM Moreno

This study evaluates the impact of LiDAR pulse density thinning on forest structure characterization using the LadderFuelsR framework. High-resolution LiDAR data were collected from Mediterranean forest sites with diverse structures and fire histories and systematically thinned to simulate varying pulse densities. Key forest metrics, including leaf area density (LAD), leaf area index (LAI), canopy base height, depth, height of fuel layers, and inter-layer distances, were derived to assess changes at different thinning levels.

# Getting Started


#SECTION 1. THINNING, CHM, SEGMENTATION, LAD PROFILES
#1.THINNING LIDARPOD FILES (CAUTION THE THINNED DATA CAN BE CHANGE IF YOU USE LIDR OR LASTOOLS AND IF USE -HIGHEST or -LOWEST parameter in LasTools)
```{r pressure, echo=FALSE}

############################
### USING LASTOOLS (in a loop for checking thinning process consistency)
############################

![] ("https://raw.githubusercontent.com/olgaviedma/Pulse-density-effects-on-LadderFuelsR-metrics/master/FIGURES_TABLES/1_THINNING_LIDARPOD.png")

```

#1. CHM FOR LAS FILES
```{r pressure, echo=FALSE}

```
#2-3. CROWNS WATERSHED AND CROWNS METRICS FOR UNTHINNED LAS FILES
```{r pressure, echo=FALSE}

```
#4. CROP ALL THINNED LAS FILES WITH TREE POLYGONS derived from BENCHMARK LAS FILE (UNTHINNED)
```{r pressure, echo=FALSE}

```

#5. LAI-LAD METRICS BY TREE and point density statistics. 6. DEPURATING TREE LAD PROFILES (\>= 5 HEIGHT VALUES)
```{r pressure, echo=FALSE}

```



#SECTION 2.LadderFuelsR METRICS AND LAI METRICS (STEPS 6-17)

#18. JOINING LADDER FUEL PROPERTIES WITH CROWN POLYGONS
```{r Joining crown polygons and ladder fuels metrics, echo=TRUE, message=FALSE, warning=FALSE}


```

#19.ALL METRICS and DIFFERENCES II: (NO THINNED - REST) 
```{r Joining crown polygons and ladder fuels metrics, echo=TRUE, message=FALSE, warning=FALSE}

## Internal tables for further processing
```

#SECTION 3. CLUSTER OF ALL TREES (across thinning levels) 
#20. PCA CLUSTER AND FREQUENCY DISTRIBUTION (FIGS.2 & 4 and FIG.S3)
```{r pressure, echo=FALSE}

###############################
## INPUT FOR TABLE S2: V.TEST (VARIABLES IMPORTANCE) AND CLUSTERS PERFORMANCE
###################################


################################################
# FIG.2: PLOT CLUSTERS IN PCA AXES 
################################################


#######################################################
## FIG.4: PLOT OF CLUSTERS DISTRIBUTION (ALL SITES TOGETHER)  
#######################################################


#######################################################
## FIG. S3: PLOT OF CLUSTERS DISTRIBUTION FACET BY SITES 
#######################################################


#######################################################
##  OTHER RESULTS DERIVED FROM CLUSTER PCA  
#######################################################
# CORRELATION OF VARIABLES WITH THE PCA DIMENSIONS


```
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

