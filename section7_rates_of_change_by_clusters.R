
#**Impact of LiDAR Pulse Density on Forest Fuels Metrics derived using LadderFuelsR**

#**Authors: Olga Viedma and JM Moreno**

#This study evaluates the impact of LiDAR pulse density thinning on forest structure characterization using the LadderFuelsR framework. 
#High-resolution LiDAR data were collected from Mediterranean forest sites with diverse structures and fire histories and systematically thinned to simulate varying pulse densities. 
#Key forest metrics, including leaf area density (LAD), leaf area index (LAI), canopy base height, depth, height of fuel layers, and inter-layer distances, were derived to assess changes at different thinning levels.
#This seventh section deals with the estimation of the rate of change [RCH] of the main fuel properties at the cluster level. 
#For doing it, we tracked the trees assigned to the clusters identified at full resolution across the different thinning levels. 
#Following, regression models were applied to each tree individually or to median values (at each thinned level) by clusters to get the median RCH from the best models.

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


###############################################################
#31.ASSIGN THE "NO_THINNED" CLUSTERS (CLUST_NT) TO THE OTHER THINNING LEVELS (FIG. 6)
###############################################################

library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra) # for grid.arrange
library(changepoint)

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


# Check the levels after reordering
levels(all_data_lidarpod$thin_level)
table(all_data_lidarpod$thin_level)

levels(all_data_lidarpod$sites)
table(all_data_lidarpod$sites)

all_data_lidarpod1 <- all_data_lidarpod


## ASSIGN THE CLUSTERS at "NO_THINNED" LEVEL"  TO THE OTHER THINNING LEVELS 

# 1. Filter data for NO_THINNED thin_level level
no_thinned_data <- all_data_lidarpod1[all_data_lidarpod1$thin_level == "NO_THINNED", ]

# 2. Group by TreeID and cluster to create a map of TreeID to cluster for NO_THINNED
library(dplyr)

# Create a mapping of TreeID to cluster for NO_THINNED
tree_cluster_map <- no_thinned_data %>%
  dplyr::select(treeID1, sites, cluster) %>%
  distinct()

# 3. Join this map with the original dataset, matching on TreeID
all_data_with_nt_cluster <- all_data_lidarpod1 %>%
  left_join(tree_cluster_map, by = c("treeID1", "sites"), suffix = c("", "_nt"))

# 4. Rename the new variable for clarity
all_data_with_nt_cluster <- all_data_with_nt_cluster %>%
  dplyr::rename(clust_nt = cluster_nt)

# The result is the original data with the added clust_nt column, where available
head(all_data_with_nt_cluster)

# Reorder factor levels
all_data_with_nt_cluster$thin_level <- factor(all_data_with_nt_cluster$thin_level, levels = new_order)
# Sort the dataframe by the thin_level column
all_data_with_nt_cluster <- all_data_with_nt_cluster[order(all_data_with_nt_cluster$thin_level), ]


write.table(all_data_with_nt_cluster, file.path(cluster_dir, "ALL_METRICS_all_sites_CLUST_PCA4_CLUST_NT.txt"), sep = "\t", row.names = FALSE)


## FIG. S3. PLOT DISTRIBUTION OF THE CLUSTERS IN EACH THINNED LEVEL ON THE BENCHMARK CLUSTER (CLUSTER_NO THINNED)

cluster_dir <- "RESULTS/CLUSTER"

all_data_with_nt_cluster <- read.table(file.path(cluster_dir, "ALL_METRICS_all_sites_CLUST_PCA4_CLUST_NT.txt"), sep = "\t", header = TRUE)

all_data_with_nt_cluster$treeID <- factor(all_data_with_nt_cluster$treeID,
                                          levels=unique(all_data_with_nt_cluster$treeID))
all_data_with_nt_cluster$thin_level <- factor(all_data_with_nt_cluster$thin_level,
                                              levels=unique(all_data_with_nt_cluster$thin_level))
all_data_with_nt_cluster$sites <- factor(all_data_with_nt_cluster$sites,
                                         levels=unique(all_data_with_nt_cluster$sites))
all_data_with_nt_cluster$cluster <- factor(all_data_with_nt_cluster$cluster,
                                           levels=unique(all_data_with_nt_cluster$cluster))
all_data_with_nt_cluster$clust_nt <- factor(all_data_with_nt_cluster$clust_nt,
                                            levels=unique(all_data_with_nt_cluster$clust_nt))

all_data_with_nt_cluster$thin_level <- factor(all_data_with_nt_cluster$thin_level, levels = new_order)
# Reorder the rows based on 'Variable'
all_data_with_nt_cluster <- all_data_with_nt_cluster[order(all_data_with_nt_cluster$thin_level), ]


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
levels(all_data_with_nt_cluster$thin_level) <- new_names

all_data_with_nt_cluster_sum <- all_data_with_nt_cluster %>%
  group_by(thin_level, clust_nt, cluster) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)

write.table (all_data_with_nt_cluster_sum, file.path(cluster_dir, "freq_distribution_cluster_PCA4_on_CLUST_NT_by_thinning_levels.txt"), sep = "\t", row.names =  F)


cluster_order<- c("1", "2", "3", "4")
all_data_with_nt_cluster_sum$clust_nt<- factor(all_data_with_nt_cluster_sum$clust_nt, levels=cluster_order )
all_data_with_nt_cluster_sum$cluster<- factor(all_data_with_nt_cluster_sum$cluster, levels=cluster_order )
all_data_with_nt_cluster_sum <- all_data_with_nt_cluster_sum[order(all_data_with_nt_cluster_sum$clust_nt), ]

names_clusters <- c("cluster 1", "cluster 2", "cluster 3", "cluster 4")

names_clusters_nt <- c("C1_NTH", "C2_NTH", "C3_NTH", "C4_NTH")

# Rename levels
levels(all_data_with_nt_cluster_sum$clust_nt) <- names_clusters_nt
levels(all_data_with_nt_cluster_sum$cluster) <- names_clusters

# Define color mapping for the clusters
cluster_colors <- c("cluster 1" = "orange", "cluster 2" = "brown", "cluster 3" = "darkgreen", "cluster 4" = "lightblue")

plot_file4 <- file.path(cluster_dir,"FIG_6_plot_distribution_clusters_by_CLUST_NT_by_thin_levels.png")

pp <- ggplot(all_data_with_nt_cluster_sum, aes(x = factor(thin_level), y = percentage, fill = factor(cluster))) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ clust_nt) +
  labs(x = "Thinning Levels", y = "Percentage", fill = "Cluster",
       title = "Percentage of Cluster Levels") +
  #coord_flip() +  # Flip coordinates for horizontal bars
  theme_bw()+
  scale_fill_manual(values = cluster_colors) +  # Define custom colors for clusters
  #theme_minimal() +  # Use a minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12,color = "black"),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 16, color = "black"),
        axis.title = element_text(size = 16,color = "black"),
        legend.title = element_text(size = 16,color = "black"),
        legend.text = element_text(size = 16,color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5,color = "black"),
        strip.text = element_text(size = 16,color = "black"))  # Increase the size of facet labels

# Display the plot
print(pp)

# Save the plot to file
ggsave(filename = plot_file4, plot = pp, width = 12, height = 8, dpi = 300)



###############################################################
#32. MEDIAN VALUES CLUSTERS (CLUST_NT)  
###############################################################

cluster_dir <- "RESULTS/CLUSTER"

all_data_with_nt_cluster <- read.table(file.path(cluster_dir, "ALL_METRICS_all_sites_CLUST_PCA4_CLUST_NT.txt"), sep = "\t", header = TRUE)

all_data_with_nt_cluster$treeID <- factor(all_data_with_nt_cluster$treeID,
                                          levels=unique(all_data_with_nt_cluster$treeID))
all_data_with_nt_cluster$thin_level <- factor(all_data_with_nt_cluster$thin_level,
                                              levels=unique(all_data_with_nt_cluster$thin_level))
all_data_with_nt_cluster$sites <- factor(all_data_with_nt_cluster$sites,
                                         levels=unique(all_data_with_nt_cluster$sites))
all_data_with_nt_cluster$cluster <- factor(all_data_with_nt_cluster$cluster,
                                           levels=unique(all_data_with_nt_cluster$cluster))
all_data_with_nt_cluster$clust_nt <- factor(all_data_with_nt_cluster$clust_nt,
                                            levels=unique(all_data_with_nt_cluster$clust_nt))

all_data_with_nt_cluster$thin_level <- factor(all_data_with_nt_cluster$thin_level, levels = new_order)
# Reorder the rows based on 'Variable'
all_data_with_nt_cluster <- all_data_with_nt_cluster[order(all_data_with_nt_cluster$thin_level), ]


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
levels(all_data_with_nt_cluster$thin_level) <- new_names


#####################################3

filtered_data_clust <- all_data_with_nt_cluster %>%
  dplyr::select(all_of(c(responses, "thin_level", "clust_nt"))) %>%
  group_by(thin_level, clust_nt) %>%
  dplyr::summarize(across(where(is.numeric),
                          list(median = ~median(.x, na.rm = TRUE),
                               sd = ~sd(.x, na.rm = TRUE),
                               lower_ci = ~quantile(.x, 0.5, na.rm = TRUE),
                               upper_ci = ~quantile(.x, 0.95, na.rm = TRUE)),
                          .names = "{.col}_{.fn}"),
                   .groups = 'drop') 



# Select only the median columns for pivoting
filtered_data_clust0 <- filtered_data_clust %>%
  dplyr::select(thin_level, clust_nt, ends_with(c("_median", "_sd")))


filtered_data_clust0_long <- filtered_data_clust0 %>%
  pivot_longer(
    cols = -c(clust_nt, thin_level),                     # All columns except 'cluster_nt'
    names_to = c("Variable", "ValueType"),   # Split the variable name and the type (median/sd)
    names_pattern = "(.*)_(median|sd)",      # Capture variable name and either median or sd
    values_to = "Value"                      # Name the new column for values
  )

filtered_data_clust0_wide <- filtered_data_clust0_long %>%
  pivot_wider(
    names_from = ValueType,  # Move 'median' and 'sd' to new columns
    values_from = Value      # Values from the 'Value' column
  )

# Check the result
filtered_data_clust0_wide

write.table (filtered_data_clust0_wide, file.path(cluster_dir, "STATISTICS_VARIABLES_BY_THINNED_LEVELS_NT_CLUSTER_WIDE.txt"), sep = "\t", row.names =  F)

################################3

filtered_data_clust_median <- filtered_data_clust %>%
  dplyr::select(thin_level, clust_nt, ends_with("_median"))

filtered_data_clust_median$clust_nt <- factor(filtered_data_clust_median$clust_nt)

write.table (filtered_data_clust_median, file.path(cluster_dir, "MEDIAN_VARIABLES_BY_THINNED_LEVELS_NT_CLUSTER_WIDE.txt"), sep = "\t", row.names =  F)


##############################3

filtered_data_clust2 <- all_data_with_nt_cluster %>%
  dplyr::select(all_of(c(responses, "clust_nt"))) %>%
  group_by(clust_nt) %>%
  dplyr::summarize(across(where(is.numeric),
                          list(median = ~median(.x, na.rm = TRUE),
                               sd = ~sd(.x, na.rm = TRUE),
                               lower_ci = ~quantile(.x, 0.05, na.rm = TRUE),
                               upper_ci = ~quantile(.x, 0.95, na.rm = TRUE)),
                          .names = "{.col}_{.fn}"),
                   .groups = 'drop') 

# Select only the median columns for pivoting
filtered_data_clust2a <- filtered_data_clust2 %>%
  dplyr::select(clust_nt, ends_with("_median"))


filtered_data_clust_long2 <- filtered_data_clust2a %>%
  pivot_longer(
    cols = -clust_nt,          # All columns except 'cluster'
    names_to = c("Variable", "ValueType"),
    names_pattern = "(.*)_(median)",  # Use regex to capture two groups
    values_to = "Value"
  )

# View the long format data
print(filtered_data_clust_long2)

# Optionally, convert 'variable' column to factor to keep the order
filtered_data_clust_long2 <- filtered_data_clust_long2 %>%
  mutate(Variable = factor(Variable, levels = unique(Variable)))

filtered_data_clust_long2<- filtered_data_clust_long2[,c(1:2,4)]

names(filtered_data_clust_long2) <- c("clust_nt",  "Variable", "mediana")

# View the long format data
print(filtered_data_clust_long2)
filtered_data_clust_long2 <- data.frame(filtered_data_clust_long2)


# Select only the median columns for pivoting
filtered_data_clust3 <- filtered_data_clust2 %>%
  dplyr::select(clust_nt, ends_with(c("_median", "_sd")))


filtered_data_clust_long3 <- filtered_data_clust3 %>%
  pivot_longer(
    cols = -clust_nt,                      # All columns except 'cluster_nt'
    names_to = c("Variable", "ValueType"),   # Split the variable name and the type (median/sd)
    names_pattern = "(.*)_(median|sd)",      # Capture variable name and either median or sd
    values_to = "Value"                      # Name the new column for values
  )

filtered_data_clust_wide <- filtered_data_clust_long3 %>%
  pivot_wider(
    names_from = ValueType,  # Move 'median' and 'sd' to new columns
    values_from = Value      # Values from the 'Value' column
  )

# Check the result
filtered_data_clust_wide

write.table (filtered_data_clust_wide, file.path(cluster_dir, "MEDIAN_VARIABLES_ONLY_NT_CLUSTER_WIDE.txt"), sep = "\t", row.names =  F)


###############################################################
#33. RATES OF CHANGE BY CLUSTERS from MEDIAN VALUES (USING THE BENCHMARK: CLUST_NT) (TABLE S6 and AND FIGS. S7-S11)
###############################################################

library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra) # for grid.arrange
library(changepoint)

cluster_dir <- "RESULTS/CLUSTER"

all_data_with_nt_cluster <- read.table(file.path(cluster_dir,"ALL_METRICS_all_sites_CLUST_PCA4_CLUST_NT.txt"), sep = "\t", header = TRUE)

# Reorder the thin_level levels to be in reverse order
all_data_with_nt_cluster$thin_level <- factor(all_data_with_nt_cluster$thin_level)
all_data_with_nt_cluster$clust_nt <- factor(all_data_with_nt_cluster$clust_nt)

table(all_data_with_nt_cluster$sites, all_data_with_nt_cluster$thin_level)

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

all_data_with_nt_cluster$thin_level <- factor(all_data_with_nt_cluster$thin_level, levels = new_order)
all_data_with_nt_cluster <- all_data_with_nt_cluster[order(all_data_with_nt_cluster$thin_level), ]


# Check the levels after reordering
levels(all_data_with_nt_cluster$thin_level)
table(all_data_with_nt_cluster$thin_level)

levels(all_data_with_nt_cluster$sites)
table(all_data_with_nt_cluster$sites)


## FITTED MODELS ON MEDIAN VARIABLES BY thin_level + CLUSTER 

library(dplyr)
library(ggplot2)
library(tidyr)

cluster_dir <- "RESULTS/CLUSTER"

filtered_data_clust_median <- read.table(file.path(cluster_dir,"MEDIAN_VARIABLES_BY_THINNED_LEVELS_NT_CLUSTER_WIDE.txt"), sep = "\t", header = TRUE)

# Check the structure of the new dataframe
str(filtered_data_clust_median)

filtered_data_clust_median$clust_nt <- factor(filtered_data_clust_median$clust_nt)


calculate_r_squared <- function(actual_values, predicted_values) {
  # Calculate the residuals
  residuals <- actual_values - predicted_values
  
  # Sum of squares of residuals
  ss_res <- sum(residuals^2)
  
  # Total sum of squares (proportional to the variance of the data)
  ss_tot <- sum((actual_values - mean(actual_values))^2)
  
  # R-squared calculation
  r_squared <- 1 - (ss_res / ss_tot)
  
  return(r_squared)
} ## not used

calculate_adjusted_r_squared <- function(actual_values, predicted_values, n_predictors) {
  # Calculate the residuals
  residuals <- actual_values - predicted_values
  
  # Sum of squares of residuals
  ss_res <- sum(residuals^2)
  
  # Total sum of squares (proportional to the variance of the data)
  ss_tot <- sum((actual_values - mean(actual_values))^2)
  
  # R-squared calculation
  r_squared <- 1 - (ss_res / ss_tot)
  
  # Number of observations
  n <- length(actual_values)
  
  # Adjusted R-squared calculation
  adjusted_r_squared <- 1 - ((1 - r_squared) * (n - 1)) / (n - n_predictors - 1)
  
  return(list(r_squared = r_squared, adjusted_r_squared = adjusted_r_squared))
}

compute_fitted_values <- function(x, intercept, slope, model_type) {
  if (model_type == "exp") {
    return(intercept * exp(slope * x))
  } else if (model_type == "log") {
    return(intercept + slope * log(x))
  } else if (model_type == "linear") {
    return(intercept + slope * x)
  } else {
    stop("Invalid model type.")
  }
}

# Define the format_formula function
format_formula <- function(intercept, slope, model_type) {
  if (model_type == "exp") {
    return(paste0("Exp Model: y = ", round(intercept, 2), " * exp(", round(slope, 2), " * x)"))
  } else if (model_type == "log") {
    return(paste0("Log Model: y = ", round(intercept, 2), " + ", round(slope, 2), " * log(x)"))
  } else if (model_type == "linear") {
    return(paste0("Linear Model: y = ", round(intercept, 2), " + ", round(slope, 2), " * x"))
  } else {
    stop("Invalid model type.")
  }
}

# Function to fit model and plot regression for multiple variables
fit_and_plot_regression <- function(original_df, variable, model_type) {
  # Ensure original_df is a data frame
  if (!is.data.frame(original_df)) {
    stop("The input original_df is not a data frame.")
  }
  
  # Ensure the variable column exists
  if (!variable %in% colnames(original_df)) {
    stop(paste("The specified variable", variable, "does not exist in the data frame."))
  }
  
  # Convert variable to numeric, if possible
  if (!is.numeric(original_df[[variable]])) {
    original_df[[variable]] <- as.numeric(original_df[[variable]])
  }
  
  # Handle NA values and check if there are enough valid points for regression
  original_df <- original_df[!is.na(original_df[[variable]]), ]
  if (nrow(original_df) < 2) {
    warning(paste("Not enough valid data points for variable", variable, "to fit a model."))
    return(NULL)
  }
  
  # Create a sequence for x values
  x <- seq_along(original_df[[variable]])
  
  # Fit the model based on the model_type parameter
  if (model_type == "linear") {
    model <- tryCatch({
      lm(original_df[[variable]] ~ x)
    }, error = function(e) {
      warning(paste("Failed to fit linear model for variable", variable))
      return(NULL)
    })
  } else if (model_type == "exp") {
    model <- tryCatch({
      nls(
        original_df[[variable]] ~ a * exp(b * x),
        start = list(a = max(original_df[[variable]], na.rm = TRUE), b = 0.1),
        algorithm = "port",  # Use a more robust algorithm
        lower = c(a = 0, b = 0),  # Define lower bounds
        upper = c(a = Inf, b = Inf)  # Define upper bounds
      )
    }, error = function(e) {
      warning(paste("Failed to fit exponential model for variable", variable))
      return(NULL)
    })
  } else if (model_type == "log") {
    model <- tryCatch({
      nls(
        original_df[[variable]] ~ a + b * log(x),
        start = list(a = min(original_df[[variable]], na.rm = TRUE), b = 1),
        algorithm = "port",  # Use a more robust algorithm
        lower = c(a = -Inf, b = -Inf),
        upper = c(a = Inf, b = Inf)
      )
    }, error = function(e) {
      warning(paste("Failed to fit logarithmic model for variable", variable))
      return(NULL)
    })
  } else {
    stop("Invalid model type.")
  }
  
  # If model fitting fails, return NULL
  if (is.null(model)) {
    return(NULL)
  }
  
  # Extract coefficients
  intercept <- coef(model)[1]
  slope <- coef(model)[2]
  
  # Create fitted values
  fitted_values <- predict(model)
  
  # Create plot
  plot <- ggplot(original_df, aes(x = x, y = .data[[variable]])) +
    geom_point(aes(color = "Original Values"), size = 2) +
    geom_line(aes(y = fitted_values, color = "Fitted Values")) +
    labs(title = paste( variable, "_", model_type, "model"),
         x = "Index",
         y = "Values",
         subtitle = paste("Intercept:", round(intercept, 2), "Slope:", round(slope, 2))) +
    scale_color_manual(values = c("Original" = "blue", "Fitted" = "red")) +
    theme_bw() +  # White background with theme_bw()
    theme(
      axis.line = element_line(color = "black"),
      axis.text.x = element_text(angle = 90, hjust = 1, size = 12, color = "black"),
      axis.text.y = element_text(size = 14, color = "black"),
      axis.title = element_text(size = 14, color = "black"),
      legend.title = element_text(size = 14, color = "black"),
      legend.text = element_text(size = 14, color = "black"),
      plot.title = element_text(size = 16, hjust = 0.5, color = "black"),
      strip.text = element_text(size = 15, color = "black")
    )
  
  
  return(plot)
}


# Initialize lists to store results across all vuelos
all_exp_results <- list()
all_log_results <- list()
all_linear_results <- list()

# Iterate through each cluster and variable in the data

for (cluster_level in unique(filtered_data_clust_median$clust_nt)) {
  
  # Filter data for the current cluster
  
  cluster_data <- filtered_data_clust_median[filtered_data_clust_median$clust_nt == cluster_level, ]
  
  # Initialize lists to store results and plots for each model type
  vuelo_exp_results <- list()
  vuelo_log_results <- list()
  vuelo_linear_results <- list()
  
  exp_plots <- list()
  log_plots <- list()
  linear_plots <- list()
  
  # Iterate over each numeric column for modeling
  for (col_name in names(cluster_data)) {
    col <- cluster_data[[col_name]]
    
    if (is.numeric(col) && length(col) == 10 && all(!is.na(col) & is.finite(col))) {
      if (any(col != 0)) {
        x <- seq_along(col)
        
        # Fitting exponential model
        exp_model <- tryCatch({
          nls(col ~ a * exp(b * x), start = list(a = max(col, na.rm = TRUE), b = 0.1))
        }, error = function(e) NULL)
        
        if (!is.null(exp_model)) {
          coefficients <- coef(exp_model)
          R_squared <- calculate_adjusted_r_squared(col, predict(exp_model), n_predictors = 1)
          vuelo_exp_results[[length(vuelo_exp_results) + 1]] <- data.frame(
            Variable = col_name, Slope = coefficients[2], Intercept = coefficients[1], 
            R_squared = R_squared$adjusted_r_squared, clust_nt = cluster_level, model = "exp"
          )
        }
        
        # Fitting logarithmic model
        log_model <- tryCatch({
          nls(col ~ a + b * log(x), start = list(a = min(col, na.rm = TRUE), b = 1))
        }, error = function(e) NULL)
        
        if (!is.null(log_model)) {
          coefficients <- coef(log_model)
          R_squared <- calculate_adjusted_r_squared(col, predict(log_model), n_predictors = 1)
          vuelo_log_results[[length(vuelo_log_results) + 1]] <- data.frame(
            Variable = col_name, Slope = coefficients[2], Intercept = coefficients[1], 
            R_squared = R_squared$adjusted_r_squared, clust_nt = cluster_level, model = "log"
          )
        }
        
        # Fitting linear model
        linear_model <- tryCatch({
          lm(col ~ x)
        }, error = function(e) NULL)
        
        if (!is.null(linear_model)) {
          coefficients <- coef(linear_model)
          R_squared <- calculate_adjusted_r_squared(col, predict(linear_model), n_predictors = 1)
          vuelo_linear_results[[length(vuelo_linear_results) + 1]] <- data.frame(
            Variable = col_name, Slope = coefficients[2], Intercept = coefficients[1], 
            R_squared = R_squared$adjusted_r_squared, clust_nt = cluster_level, model = "lin"
          )
        }
      }
    }
  }
  
  # Append model results for each variable
  all_exp_results <- rbind(all_exp_results, do.call(rbind, vuelo_exp_results))
  all_linear_results <- rbind(all_linear_results, do.call(rbind, vuelo_linear_results))
  all_log_results <- rbind(all_log_results, do.call(rbind, vuelo_log_results))
  
  # Plot after results collection for each variable in this cluster
  for (model_type in c("exp", "log", "linear")) {
    model_results <- get(paste0("vuelo_", model_type, "_results"))
    
    if (length(model_results) > 0) {
      # Iterate over the variables in model_results for plotting
      for (result in model_results) {
        var <- result$Variable
        plot <- fit_and_plot_regression(cluster_data, var, model_type)
        
        # Append to corresponding plot list
        if (model_type == "exp") {
          exp_plots[[length(exp_plots) + 1]] <- plot
        } else if (model_type == "log") {
          log_plots[[length(log_plots) + 1]] <- plot
        } else if (model_type == "linear") {
          linear_plots[[length(linear_plots) + 1]] <- plot
        }
      }
    }
  }
  
  # Save the plots for each model type in chunks of 12
  for (model_type in c("exp", "log", "linear")) {
    plots <- get(paste0(model_type, "_plots"))
    
    if (length(plots) > 0) {
      plot_chunks <- split(plots, ceiling(seq_along(plots) / 12))
      
      for (i in seq_along(plot_chunks)) {
        combined_plot <- cowplot::plot_grid(plotlist = plot_chunks[[i]], ncol = 3)
        ggsave(
          filename = file.path(cluster_dir, paste0("regress_plot_NT_cluster_", cluster_level, "_model_", model_type, "_chunk_", i, ".png")),
          plot = combined_plot,
          width = 14, height = 8
        )
      }
    }
  }
  
}

# Combine results across all vuelos
combined_exp_results <- do.call(rbind, all_exp_results)
combined_log_results <- do.call(rbind, all_log_results)
combined_linear_results <- do.call(rbind, all_linear_results)

combined_exp_results_df <- as.data.frame(t(combined_exp_results))
combined_log_results_df <- as.data.frame(t(combined_log_results))
combined_lin_results_df <- as.data.frame(t(combined_linear_results))

# Write the combined results to a single file for each model
write.table(combined_lin_results_df, file.path(cluster_dir, "LINEAR_FUNCTIONS_MEDIAN_VARIABLES_NT_CLUSTER.txt"), sep = "\t", row.names = FALSE)
write.table(combined_exp_results_df, file.path(cluster_dir, "EXP_FUNCTIONS_MEDIAN_VARIABLES_NT_CLUSTER.txt"), sep = "\t", row.names = FALSE)
write.table(combined_log_results_df, file.path(cluster_dir, "LOG_FUNCTIONS_MEDIAN_VARIABLES_NT_CLUSTER.txt"), sep = "\t", row.names = FALSE)




## JOINING ALL FITTED MODELs OVER MEDIAN VARIABLES (BY THINNING LEVELS + CLUSTER NT)

cluster_dir <- "RESULTS/CLUSTER"

models_list <- list.files(cluster_dir, pattern = glob2rx("*_FUNCTIONS_MEDIAN_VARIABLES_NT_CLUSTER.txt"), full.names = TRUE)
models_file <- lapply(models_list, function(X) read.table(X, sep = "\t", header = TRUE))


all_models <- do.call (rbind, models_file)
# Filter the data to exclude rows where Variable is "clust_nt"
all_models <- all_models[all_models$Variable != "clust_nt", ]

write.table(all_models, file.path(cluster_dir, "ALL_FITTED_FUNCTIONS_MEDIAN_VARIABLES_NT_CLUSTER.txt"), sep = "\t", row.names = FALSE)


## BEST FITTED MODELS BY CLUST_NT

# Find the best model for each variable

best_models <- all_models %>%
  group_by(Variable,clust_nt) %>%
  dplyr::filter(R_squared == max(R_squared)) %>%
  arrange(Variable,clust_nt)

best_models_clean <- best_models %>%
  filter(R_squared != -Inf)

write.table(best_models_clean, file.path(cluster_dir, "BEST_R2_FITTED_FUNCTIONS_MEDIAN_VARIABLES_NT_CLUSTER.txt"), sep = "\t", row.names = FALSE)


#  CALCULATE RATES OF CHANGE (SLOPE) OF BEST "MEDIAN" MODELS BY thin_level + CLUSTER

library(dplyr)

# Function to calculate percentage change for exponential models
calculate_percentage_change_exp <- function(Slope) {
  return((exp(Slope) - 1) * 100)
}

# Function to calculate the rate of change for logarithmic models
calculate_rate_of_change_log <- function(Slope) {
  return(Slope * 0.01)  # Slope gives the percentage change in response to a 1% change in the predictor
}

# Updated best_models_clean1 calculation with rate of change for different models
best_models_clean1 <- best_models_clean %>%
  mutate(rate_of_change = case_when(
    model == "lin" ~ Slope * 100,  # Linear model rate of change
    model == "exp" ~ calculate_percentage_change_exp(Slope),  # Exponential model rate of change
    model == "log" ~ calculate_rate_of_change_log(Slope),  # Logarithmic model rate of change
    TRUE ~ NA_real_  # Default to NA for other models
  )) %>%
  ungroup()  # Ungroup after the calculation

# Print the updated table with percentage change
print(best_models_clean1)

# Round numeric variables to 2 decimal places
best_models_clean1 <- best_models_clean1 %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

best_models_clean1$clust_nt <- factor(best_models_clean1$clust_nt)

best_models_clean1 <- best_models_clean1 %>% mutate(Variable = str_remove(Variable, "_median$"))  # Use regex to remove the last underscore

# Perform the left join

filtered_data_clust2a_wide <- read.table(file.path(cluster_dir,"MEDIAN_VARIABLES_ONLY_NT_CLUSTER_WIDE.txt"), sep = "\t", header = TRUE)

cluster_order<- c("1", "2", "3", "4")
best_models_clean1$clust_nt<- factor(best_models_clean1$clust_nt, levels=cluster_order )
filtered_data_clust2a_wide$clust_nt <- factor(filtered_data_clust2a_wide$clust_nt, levels = cluster_order)

best_models_clean2 <- best_models_clean1 %>%
  left_join(filtered_data_clust2a_wide, by = c("Variable", "clust_nt"))


write.table(best_models_clean2, file.path(cluster_dir, "RATES_CHANGES_BEST_R2_MEDIAN_VALUES_CLUST_NT.txt"), sep = "\t", row.names = FALSE)


# TABLE_S6: CONVERT BEST MODELS MEDIANS (BY thin_level + CLUSTER) AND RATES OF CHANGE (SLOPE) INTO A .DOCX TABLE


metrics_crowns_dir<- "CHMs/WATER_LAI/CLUSTER_ALL_DIFF_THINNED"

best_models_clean2 <- read.table(file.path(cluster_dir, "RATES_CHANGES_BEST_R2_MEDIAN_VALUES_CLUST_NT.txt"), sep = "\t", header = TRUE)

best_models_clean2$Variable <- factor(best_models_clean2$Variable,
                                      levels=unique(best_models_clean2$Variable))
best_models_clean2$model <- factor(best_models_clean2$model,
                                   levels=unique(best_models_clean2$model))
best_models_clean2$clust_nt <- factor(best_models_clean2$clust_nt,
                                      levels=unique(best_models_clean2$clust_nt))


responses <-c("F1_cbh", "F1_depth", "F1_LAD","F1_Hdepth","MXLAD_cbh", "BP_cbh", "BR_cbh",
              "MXLAD_depth", "BP_depth", "MXLAD_dist", "BP_dist",
              "MXLAD_Hdepth","BP_Hdepth", "MAX_height", "total_LAI", "understory_LAI")


best_models_clean2 <- best_models_clean2 %>%
  dplyr::filter(Variable %in% responses)

# Reorder factor levels for the Variable column based on the desired order
best_models_clean2$Variable <- factor(best_models_clean2$Variable, levels = responses)

# Reorder the rows based on 'Variable'
best_models_clean2 <- best_models_clean2[order(best_models_clean2$Variable), ]


library(officer)
library(flextable)

best_models_clean2_formatted <- best_models_clean2 %>%
  dplyr::mutate(
    across(
      .cols = where(is.numeric) & !all_of("R_squared"), 
      .fns = ~ as.numeric(formatC(., format = "f", digits = 1))
    ),
    R_squared = as.numeric(formatC(R_squared, format = "f", digits = 2))
  )


best_models_clean3 <- best_models_clean2_formatted %>%
  # Combine 'median' and 'sd' into the 'mediana' column with formatted values
  mutate(median = paste0(median, " ± ", sd)) %>%
  arrange(Variable) %>%
  dplyr::select(Variable, clust_nt, model, Intercept, Slope, R_squared, rate_of_change, median)


# Create a Word document
doc <- read_docx()

# Convert the tibble to a flextable
ft <- flextable(best_models_clean3)
autofit(ft)  # Automatically adjusts the column width to fit the content

# Add the flextable to the document
doc <- doc %>% 
  body_add_flextable(ft)

# Save the document
print(doc, target = file.path(cluster_dir,"TABLE_S6_table_best_models_RATE_CHANGE_MEDIAN_VALUES_CLUSTER_NT.docx"))



# FIGS. S7-S11.PLOT FITTED MODELS AND RATES OF CHANGES OVER MEDIAN VALUES USING CLUSTER_NT

output_dir<- "CHMs/WATER_LAI/CLUSTER_ALL_DIFF_THINNED"

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(ggrepel)

best_models_clean2 <- read.table(file.path(cluster_dir,"RATES_CHANGES_BEST_R2_MEDIAN_VALUES_CLUST_NT.txt"), sep = "\t", header = TRUE)

best_models_clean2$Variable <-factor(best_models_clean2$Variable )
best_models_clean2$clust_nt <-factor(best_models_clean2$clust_nt )

# Check the updated levels
levels(best_models_clean2$Variable)

responses <-c("F1_cbh", "F1_depth", "F1_LAD","F1_Hdepth","MXLAD_cbh", "BP_cbh", "BR_cbh",
              "MXLAD_depth", "BP_depth", "MXLAD_dist", "BP_dist",
              "MXLAD_Hdepth","BP_Hdepth", "MAX_height", "total_LAI", "understory_LAI")

# Reorder factor levels
best_models_clean2$Variable <- factor(best_models_clean2$Variable, levels = responses)
# Reorder the rows based on 'Variable'
best_models_clean2 <- best_models_clean2[order(best_models_clean2$Variable), ]

##############3

all_data_with_nt_cluster <- read.table(file.path(cluster_dir,"ALL_METRICS_all_sites_CLUST_PCA4_CLUST_NT.txt"), sep = "\t", header = TRUE)

# Reorder the thin_level levels to be in reverse order
all_data_with_nt_cluster$thin_level <- factor(all_data_with_nt_cluster$thin_level)
all_data_with_nt_cluster$clust_nt <- factor(all_data_with_nt_cluster$clust_nt)

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
all_data_with_nt_cluster$thin_level <- factor(all_data_with_nt_cluster$thin_level, levels = new_order)
all_data_with_nt_cluster <- all_data_with_nt_cluster[order(all_data_with_nt_cluster$thin_level), ]

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
levels(all_data_with_nt_cluster$thin_level) <- new_names

# Define your variable groups
variable_groups <- list(
  underst =c("F1_cbh", "F1_depth", "F1_LAD", "F1_Hdepth"),
  Hcbh = c( "MXLAD_cbh", "BP_cbh", "BR_cbh"),
  Depth = c("MXLAD_depth", "BP_depth","MXLAD_dist", "BP_dist"),
  Hdepth = c("MXLAD_Hdepth", "BP_Hdepth"),
  Other = c("MAX_height", "total_LAI", "understory_LAI")
)


# Define the desired order of clusters
cluster_order <- c("1", "2", "3", "4")

# Initialize an empty list to store combined plots for each variable group
combined_plots_group <- list()

# Loop over each variable group
for (group_name in names(variable_groups)) {
  
  # Get variables for the current group
  vars_to_plot <- variable_groups[[group_name]]
  
  # Initialize a list to store individual plots for each variable
  individual_plots <- list()
  
  # Loop over each variable in the group
  for (var in vars_to_plot) {
    
    # Filter and summarize data for the current variable
    filtered_data <- all_data_with_nt_cluster %>%
      dplyr::select(all_of(c(var, "thin_level", "clust_nt"))) %>%
      group_by(thin_level, clust_nt) %>%
      summarize(
        media = median(.data[[var]], na.rm = TRUE),
        sd = sd(.data[[var]], na.rm = TRUE),
        lower_ci = quantile(.data[[var]], 0.25, na.rm = TRUE),
        upper_ci = quantile(.data[[var]], 0.75, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      mutate(
        lower_bound = lower_ci,
        upper_bound = upper_ci
      )
    
    if (nrow(filtered_data) > 0) {
      # Order clusters
      filtered_data$clust_nt <- factor(filtered_data$clust_nt, levels = cluster_order)
      
      # Initialize annotations
      annotations <- data.frame(thin_level = character(), value_at_first_thin_level = numeric(),
                                R_squared = numeric(), Variable = character(),
                                clust_nt = character(), stringsAsFactors = FALSE)
      
      # Extract R-squared and rate of change for annotations
      for (cl in levels(filtered_data$clust_nt)) {
        filtered_data1 <- filtered_data %>% filter(clust_nt == cl)
        rsq_values <- best_models_clean2 %>%
          filter(Variable == var, clust_nt == cl) %>%
          distinct() %>%
          mutate(R_squared = as.numeric(R_squared), rate_of_change = as.numeric(rate_of_change))
        
        if (nrow(rsq_values) > 0) {
          first_value <- filtered_data1 %>%
            summarize(first_value = max(media[!is.na(media)]),
                      first_thin_level = first(thin_level), .groups = 'drop')
          annotations <- rbind(annotations, data.frame(
            thin_level = first_value$first_thin_level,
            value_at_first_thin_level = first_value$first_value,
            R_squared = round(rsq_values$R_squared[1], 2),
            rate_of_change = round(rsq_values$rate_of_change[1], 2),
            Variable = var, clust_nt = cl
          ))
        }
      }
      
      # Plot with adjustments
      plot_variable <- ggplot(filtered_data, aes(x = thin_level, y = media, group = clust_nt, color = clust_nt)) +
        geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound, fill = clust_nt), alpha = 0.3) +
        geom_line(size = 1) +
        geom_point(size = 2) +
        labs(title = var, x = "Thinning levels", y = "Median ± P25-75") +
        scale_color_manual(values = c("1" = "orange", "2" = "brown", "3" = "darkgreen", "4" = "lightblue")) +
        scale_fill_manual(values = c("1" = "orange", "2" = "brown", "3" = "darkgreen", "4" = "lightblue")) +
        theme_bw() +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1, size = 14, color = "black"),
          axis.text.y = element_text(size = 14, color = "black"),
          axis.title = element_text(size = 14, color = "black"),
          plot.title = element_text(size = 16, hjust = 0.5, color = "black"),
          strip.text = element_text(size = 14, color = "black")
        ) +
        facet_wrap(~ clust_nt, ncol = 2, scales = "free_y")
      
      # Add annotations
      if (nrow(annotations) > 0) {
        annotations <- annotations %>%
          mutate(annotation_y = value_at_first_thin_level + 0.2 * (max(filtered_data$media, na.rm = TRUE) - min(filtered_data$media, na.rm = TRUE)))
        
        plot_variable <- plot_variable +
          geom_text_repel(
            data = annotations,
            aes(x = thin_level, y = annotation_y, label = paste0("R²: ", R_squared, " RCH: ", rate_of_change)),
            size = 5, color = "black", box.padding = 0.4, point.padding = 0.4, segment.color = "grey50"
          )
      }
      
      # Store the plot
      individual_plots[[var]] <- plot_variable
    }
  }
  
  # Pad empty plots for groups with <4 variables
  while (length(individual_plots) < 4) {
    individual_plots[[paste0("empty_", length(individual_plots) + 1)]] <- ggplot() +
      theme_void() +
      labs(title = "", x = NULL, y = NULL)
  }
  
  # Combine plots for the group
  combined_plot <- wrap_plots(individual_plots, ncol = 2) +
    plot_annotation(title = paste("Group:", group_name))
  
  # Save the combined plot
  ggsave(
    filename = file.path(cluster_dir, paste0("plot_RCH_quantiles_2575_CLUSTER_NT_", group_name, ".png")), 
    plot = combined_plot, 
    width = 12, height = 9, dpi = 300
  )
  
  # Store in combined_plots_group
  combined_plots_group[[group_name]] <- combined_plot
}


###############################################################
#34.COMPARISON RATES OF CHANGE OF INDIVIDUAL TREES vs. MEDIAN VALUES (CLUST_NT) (TABLE S7 and FIG. S6)
###############################################################

crowns_dir <- "RESULTS/MODELS"

best_models<- read.table (file.path(crowns_dir, "BEST_R2_FITTED_ALL_FUNCTIONS_ALL_TREES.txt"), sep = "\t", header = TRUE)

best_models$treeID <- factor(best_models$treeID,levels=unique(best_models$treeID))
best_models$sites <- factor(best_models$sites,levels=unique(best_models$sites))

# Create a new variable "treeID1" by concatenating "treeID" and "sites"
best_models$treeID_site <- factor(
  paste(best_models$treeID, best_models$sites, sep = "_"),
  levels = unique(paste(best_models$treeID, best_models$sites, sep = "_"))
)

###################3

cluster_dir <- "RESULTS/CLUSTER"

all_data_lidarpod <- read.table(file.path(cluster_dir, "ALL_METRICS_all_sites_CLUST_PCA4_CLUST_NT.txt"), sep = "\t", header = TRUE)

all_data_lidarpod$treeID <- factor(all_data_lidarpod$treeID,levels=unique(all_data_lidarpod$treeID))
all_data_lidarpod$thin_level <- factor(all_data_lidarpod$thin_level,levels=unique(all_data_lidarpod$thin_level))
all_data_lidarpod$sites <- factor(all_data_lidarpod$sites,levels=unique(all_data_lidarpod$sites))
all_data_lidarpod$cluster <- factor(all_data_lidarpod$cluster,levels=unique(all_data_lidarpod$cluster))
all_data_lidarpod$clust_nt <- factor(all_data_lidarpod$clust_nt,levels=unique(all_data_lidarpod$clust_nt))

# Create a new variable "treeID1" by concatenating "treeID" and "sites"
all_data_lidarpod$treeID_site <- factor(
  paste(all_data_lidarpod$treeID1, all_data_lidarpod$sites, sep = "_"),
  levels = unique(paste(all_data_lidarpod$treeID1, all_data_lidarpod$sites, sep = "_"))
)

all_data_lidarpod_set <-all_data_lidarpod %>%
  dplyr::filter(thin_level== "NO_THINNED")

# Count occurrences of each treeID and filter for those with more than 1 row
treeID_multiple_rows <- all_data_lidarpod_set %>%
  count(treeID_site) %>%
  filter(n > 1)

# Check if there are any treeIDs with more than 1 row
if (nrow(treeID_multiple_rows) > 0) {
  print("TreeIDs with more than 1 row:")
  print(treeID_multiple_rows)
} else {
  print("No treeIDs with more than 1 row.")
}

all_models_cluster <- left_join(best_models, all_data_lidarpod_set, by = "treeID_site")
names(all_models_cluster) <- sub("\\.x$", "", names(all_models_cluster))
cols_to_drop <- c("treeID.y", "sites.y")
all_models_cluster <- all_models_cluster[, !names(all_models_cluster) %in% cols_to_drop]

all_models_cluster <- all_models_cluster %>%
  # Select the relevant columns, keeping only one 'sites' column
  dplyr::select(treeID, treeID_site, sites, Variable, Slope, Intercept, R_squared, model, clust_nt)

# Check the modified dataframe
head(all_models_cluster)

write.table (all_models_cluster, file.path(cluster_dir, "BEST_R2_FITTED_FUNCTIONS_ALL_TREES_NT_CLUSTER.txt"), sep = "\t", row.names =  F)

##############################3
# Define the desired order of responses
responses <-c("F1_cbh", "F1_depth", "F1_LAD","F1_Hdepth","MXLAD_cbh", "BP_cbh", "BR_cbh",
              "MXLAD_depth", "BP_depth", "MXLAD_dist", "BP_dist",
              "MXLAD_Hdepth","BP_Hdepth", "MAX_height", "total_LAI", "understory_LAI")

filtered_data <- all_data_lidarpod %>%
  dplyr::select(all_of(c(responses, "thin_level", "clust_nt"))) %>%
  group_by(thin_level, clust_nt) %>%
  dplyr::summarize(across(where(is.numeric),
                          list(median = ~median(.x, na.rm = TRUE),
                               sd = ~sd(.x, na.rm = TRUE),
                               lower_ci = ~median(.x, na.rm = TRUE) - 1.96 * sd(.x, na.rm = TRUE) / sqrt(n()),
                               upper_ci = ~median(.x, na.rm = TRUE) + 1.96 * sd(.x, na.rm = TRUE) / sqrt(n())),
                          .names = "{.col}_{.fn}"),
                   .groups = 'drop') 


filtered_data1 <- filtered_data %>%
  dplyr::select(thin_level, clust_nt, ends_with("_median"))



###############3

reference_values_long1 <- filtered_data1 %>%
  pivot_longer(
    cols = F1_cbh_median :understory_LAI_median,  # Specify the range of columns to pivot
    names_to = "Variable", # The new column name that will hold the variable names
    values_to = "median"    # The new column name that will hold the values
  )

reference_values_long1$Variable<-factor(reference_values_long1$Variable)
reference_values_long1$Variable <- factor(gsub("_median$", "", levels(reference_values_long1$Variable)[reference_values_long1$Variable]))

write.table (reference_values_long1, file.path(cluster_dir, "MEDIAN_ALL_TREES_NT_CLUSTER_all_THINNED_LEVELS.txt"), sep = "\t", row.names =  F)



# INPUT DATA FOR RATES OF CHANGE USING ALL TREES" BY thin_level + CLUST_NT

cluster_dir <- "RESULTS/CLUSTER"

best_models_medians<- read.table (file.path(cluster_dir, "RATES_CHANGES_BEST_R2_MEDIAN_VALUES_CLUST_NT.txt"), sep = "\t", header = TRUE)
best_models_medians$Variable<-factor(best_models_medians$Variable)
best_models_medians$model<-factor(best_models_medians$model)
best_models_medians$clust_nt<-factor(best_models_medians$clust_nt)
best_models_medians$Variable <- factor(gsub("_median$", "", levels(best_models_medians$Variable)[best_models_medians$Variable]))

# Add "_median" to all columns except "Variable" and "model"
colnames(best_models_medians) <- ifelse(colnames(best_models_medians) %in% c("Variable", "model" ,"clust_nt"),
                                        colnames(best_models_medians),
                                        paste0(colnames(best_models_medians), "_median"))

# Check the updated column names
names(best_models_medians)

best_models<- read.table (file.path(cluster_dir, "BEST_R2_FITTED_FUNCTIONS_ALL_TREES_NT_CLUSTER.txt"), sep = "\t", header = TRUE)
best_models$treeID<-factor(best_models$treeID)
best_models$treeID_site<-factor(best_models$treeID_site)
best_models$sites<-factor(best_models$sites)
best_models$Variable<-factor(best_models$Variable)
best_models$model<-factor(best_models$model)
best_models$clust_nt<-factor(best_models$clust_nt)


# Perform a join to match best_models_clean with best_models_medians on Variable and model
best_models_clean1 <- best_models %>%
  inner_join(best_models_medians,  by = c("Variable", "model", "clust_nt"))

# View the summary of the matched models
summary(best_models_clean1)

best_models_clean1$model <- factor(best_models_clean1$model)


responses <-c("F1_cbh", "F1_depth", "F1_LAD","F1_Hdepth","MXLAD_cbh", "BP_cbh", "BR_cbh",
              "MXLAD_depth", "BP_depth", "MXLAD_dist", "BP_dist",
              "MXLAD_Hdepth","BP_Hdepth", "MAX_height", "total_LAI", "understory_LAI")

best_models_clean2 <- best_models_clean1 %>%
  dplyr::filter(Variable %in% responses)


#  CALCULATE RATES OF CHANGE (SLOPE) OF BEST MODELS FROM " ALL TREES" (BY thin_level + CLUST_NT)

library(dplyr)

# Function to calculate percentage change for exponential models
calculate_percentage_change_exp <- function(Slope) {
  return((exp(Slope) - 1) * 100)
}

# Function to calculate the rate of change for logarithmic models
calculate_rate_of_change_log <- function(Slope) {
  return(Slope * 0.01)  # Slope gives the percentage change in response to a 1% change in the predictor
}

# Updated best_models_clean1 calculation with rate of change for different models
best_models_cluster1 <- best_models_clean2 %>%
  mutate(rate_of_change = case_when(
    model== "lin" ~ Slope * 100,  # Linear model rate of change
    model == "exp" ~ calculate_percentage_change_exp(Slope),  # Exponential model rate of change
    model == "log" ~ calculate_rate_of_change_log(Slope),  # Logarithmic model rate of change
    TRUE ~ NA_real_  # Default to NA for other models
  )) %>%
  ungroup()  # Ungroup after the calculation

# Print the updated table with percentage change
print(best_models_cluster1)

# Round numeric variables to 2 decimal places
best_models_cluster1 <- best_models_cluster1 %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

write.table(best_models_cluster1, file.path(cluster_dir, "RATES_CHANGES_CLUST_NT_ALL_TREES_plus_MEDIAN_MODELS.txt"), sep = "\t", row.names = FALSE)



# TABLE S7. CONVERT BEST MODELS FROM ALL TREES (BY thin_level + CLUSTER ) AND RATES OF CHANGE (SLOPE) INTO A .DOCX TABLE

library(officer)
library(flextable)


best_models_cluster2<- read.table (file.path(cluster_dir, "RATES_CHANGES_CLUST_NT_ALL_TREES_plus_MEDIAN_MODELS.txt"), sep = "\t", header = TRUE)

best_models_cluster2$treeID<-factor(best_models_cluster2$treeID)
best_models_cluster2$sites<-factor(best_models_cluster2$sites)
best_models_cluster2$Variable<-factor(best_models_cluster2$Variable)
best_models_cluster2$model<-factor(best_models_cluster2$model)
best_models_cluster2$clust_nt<-factor(best_models_cluster2$clust_nt)

#################### USING QUANTILES

summary_data3 <- best_models_cluster2 %>%
  dplyr::group_by(Variable, clust_nt) %>%
  dplyr::summarize(
    rate_change_trees1 = median(rate_of_change, na.rm = TRUE),
    Slope_trees1 = median(Slope, na.rm = TRUE),
    sd_slope_trees1 = sd(Slope, na.rm = TRUE),
    Intercept_trees1 = median(Intercept, na.rm = TRUE),
    sd_intercept_trees1 = sd(Intercept, na.rm = TRUE),
    R2_trees1 = median(R_squared, na.rm = TRUE),
    sd_R2_trees1 = sd(R_squared, na.rm = TRUE),
    R2_25th_percentile = quantile(R_squared, 0.25, na.rm = TRUE),
    R2_75th_percentile = quantile(R_squared, 0.75, na.rm = TRUE),
    mediana1 = median( median_median, na.rm = TRUE),
    sd_mediana1 = median(sd_median, na.rm = TRUE),
    
    # Calculate Q1 and Q3 for rate_of_change
    Q1_Rtch = quantile(rate_of_change, 0.05, na.rm = TRUE),
    Q3_Rtch= quantile(rate_of_change, 0.95, na.rm = TRUE)
  ) %>%
  mutate(
    # Calculate bounds using Q1 and Q3
    lower_bound = Q1_Rtch,
    upper_bound = Q3_Rtch,
    Variable_numeric = as.numeric(factor(Variable))  # Create a numeric version of Variable
  )


summary_data3<- summary_data3 %>%
  mutate(across(where(is.character), ~ as.numeric(as.character(.))))

summary_data3 <- summary_data3 %>%
  mutate(across(
    where(is.numeric) & !all_of(c("R2_trees1")),
    ~ formatC(., format = "f", digits = 1)
  )) %>%
  mutate(across(
    all_of(c("R2_trees1")),
    ~ formatC(., format = "f", digits = 2)
  ))

summary_data4 <- summary_data3 %>%
  dplyr::mutate(Intercept = paste0(Intercept_trees1, " ± ", sd_intercept_trees1),
                Slope = paste0(Slope_trees1, " ± ", sd_slope_trees1),
                R2 = R2_trees1,
                bounds_R2 = paste0("(", R2_25th_percentile, ",", R2_75th_percentile, ")"),
                #rate_change_trees = paste0(rate_change_trees1, " ± ", sd_change_trees1),
                RCH = paste0(rate_change_trees1),
                bounds_RCH = paste0("(", lower_bound, ",", upper_bound, ")"),
                median = paste0(mediana1, " ± ", sd_mediana1)) %>%  # closed parentheses correctly
  arrange(Variable) %>%
  dplyr::select(Variable, clust_nt, Intercept, Slope,  R2, bounds_R2, RCH,bounds_RCH, median)

responses <-c("F1_cbh", "F1_depth", "F1_LAD","F1_Hdepth","MXLAD_cbh", "BP_cbh", "BR_cbh",
              "MXLAD_depth", "BP_depth", "MXLAD_dist", "BP_dist",
              "MXLAD_Hdepth","BP_Hdepth", "MAX_height", "total_LAI", "understory_LAI")


summary_data4$Variable<- factor(summary_data4$Variable, levels=responses)

# Reorder the rows based on 'Variable'
summary_data4 <- summary_data4[order(summary_data4$Variable), ]

summary_data4 <- summary_data4 %>%
  mutate(across(c(R2), 
                ~ formatC(., format = "f", digits = 2)))

# Create a Word document
doc2a <- read_docx()

# Create the flextable
ft2a <- flextable(summary_data4)

# Try autofit after confirming the table is generated correctly
ft2a <- autofit(ft2a)

# Add the flextable to the document
doc2a <- doc2a %>% 
  body_add_flextable(ft2a)

# Save the document
print(doc2a, target = file.path(crowns_dir,"TABLE_S7_best_models_RATE_CHANGE_ALL_TREES_CLUSTER_NT_p2575.docx"))


## KRUSKAL-WALLIS NT_CLUSTER (R2 AND RATE OF CHANGE)

library(phia)
library(dunn.test)
library(ggplot2)
library(dplyr)
library(tidyr)
library(FSA)  # for kwAllPairsDunnTest
library(lsr)  # for etaSquared
library(gridExtra)
library(PMCMRplus)
library(multcompView)  # For multcompLetters function


all_models_cluster <- read.table( file.path(cluster_dir, "RATES_CHANGES_CLUST_NT_ALL_TREES_plus_MEDIAN_MODELS.txt"), sep = "\t", header = TRUE)

all_models_cluster$treeID <- factor(all_models_cluster$treeID,
                                    levels=unique(all_models_cluster$treeID))
all_models_cluster$treeID_site <- factor(all_models_cluster$treeID_site,
                                         levels=unique(all_models_cluster$treeID_site))
all_models_cluster$sites <- factor(all_models_cluster$sites,
                                   levels=unique(all_models_cluster$sites))
all_models_cluster$clust_nt <- factor(all_models_cluster$clust_nt,
                                      levels=unique(all_models_cluster$clust_nt))
all_models_cluster$model <- factor(all_models_cluster$model,
                                   levels=unique(all_models_cluster$model))

responses <-c("F1_cbh", "F1_depth", "F1_LAD","F1_Hdepth","MXLAD_cbh", "BP_cbh", "BR_cbh",
              "MXLAD_depth", "BP_depth", "MXLAD_dist", "BP_dist",
              "MXLAD_Hdepth","BP_Hdepth", "MAX_height", "total_LAI", "understory_LAI")


all_models_cluster <- all_models_cluster %>%
  dplyr::filter(Variable %in% responses)

# Reorder factor levels for the Variable column based on the desired order
all_models_cluster$Variable <- factor(all_models_cluster$Variable, levels = responses)

# Reorder the rows based on 'Variable'
all_models_cluster <- all_models_cluster[order(all_models_cluster$Variable), ]

unique_clusters <- unique(all_models_cluster$clust_nt)


# Define the flexible function to perform Kruskal-Wallis tests and related calculations
perform_kruskal_analysis <- function(data, variable, response_var) {
  results <- list()
  eta_results <- list()
  letters_results <- list()
  
  # Filter the data for the current variable
  data_subset <- data %>%
    dplyr::filter(variable == Variable)
  
  # Check if the response variable has all zero values
  if (all(data_subset[[response_var]] == 0)) {
    cat("Skipping", variable, "as all", response_var, "values are zero.\n")
    return(list(results = NULL, eta_results = NULL, letters_results = NULL))
  }
  
  # Check if all observations belong to a single group
  if (length(unique(data_subset$clust_nt)) == 1) {
    cat("Skipping", variable, "as all observations belong to a single group for", response_var, "\n")
    return(list(results = NULL, eta_results = NULL, letters_results = NULL))
  }
  
  
  # Perform Kruskal-Wallis test
  kw_test <- kwAllPairsDunnTest(as.formula(paste(response_var, "~ clust_nt")), data = data_subset)
  p_values <- kw_test$p.value
  valid_comparisons <- which(!is.na(as.vector(p_values)))
  
  # Generate names for valid comparisons
  named_p_values <- as.vector(p_values)[valid_comparisons]
  comparison_names <- apply(expand.grid(rownames(p_values), colnames(p_values)), 1, function(x) paste(sort(x), collapse = "-"))
  comparison_names <- comparison_names[valid_comparisons]
  names(named_p_values) <- comparison_names
  
  # Generate significance letters
  sig_letters <- multcompLetters(p.adjust(named_p_values, method = "bonferroni"))
  letter_df <- data.frame(cluster = names(sig_letters$Letters), 
                          letter = sig_letters$Letters, 
                          variable = variable, 
                          measure = response_var)
  letters_results[[variable]] <- letter_df
  
  # Transform p-values to a long format for further processing
  p_values_df <- as.data.frame(p_values)
  p_values_df$cluster <- rownames(p_values_df)
  p_values_long <- tidyr::gather(p_values_df, comparison, p_value, -cluster, na.rm = TRUE)
  p_values_long$variable <- variable
  p_values_long$measure <- response_var
  
  results[[variable]] <- p_values_long
  
  # Perform ANOVA and calculate eta squared
  aov_model <- aov(as.formula(paste(response_var, "~ clust_nt")), data = data_subset)
  eta1 <- as.data.frame(etaSquared(aov_model, type = 2))
  eta1$variable <- variable
  eta1$measure <- response_var
  
  eta_results[[variable]] <- eta1
  
  return(list(results = results, eta_results = eta_results, letters_results = letters_results))
}

# Define the levels of Variable
variable_levels <-c("F1_cbh", "F1_depth", "F1_LAD","F1_Hdepth","MXLAD_cbh", "BP_cbh", "BR_cbh",
                    "MXLAD_depth", "BP_depth", "MXLAD_dist", "BP_dist",
                    "MXLAD_Hdepth","BP_Hdepth", "MAX_height", "total_LAI", "understory_LAI")

# Specify the response variables you want to analyze
response_vars <- c("R_squared", "rate_of_change")  # Add more as needed

# Initialize lists to hold combined results
all_results_kruskal <- list()
all_results_eta <- list()
all_results_letters <- list()

# Loop through each specified variable and each response variable
for (var in variable_levels) {
  for (response_var in response_vars) {
    cat("Processing variable:", var, "for response:", response_var, "\n")
    analysis_results <- perform_kruskal_analysis(all_models_cluster, var, response_var)
    
    # Combine results if they are not NULL
    if (!is.null(analysis_results$results)) {
      if (is.null(all_results_kruskal[[var]])) {
        all_results_kruskal[[var]] <- analysis_results$results[[var]]
      } else {
        all_results_kruskal[[var]] <- rbind(all_results_kruskal[[var]], analysis_results$results[[var]])
      }
    }
    if (!is.null(analysis_results$eta_results)) {
      if (is.null(all_results_eta[[var]])) {
        all_results_eta[[var]] <- analysis_results$eta_results[[var]]
      } else {
        all_results_eta[[var]] <- rbind(all_results_eta[[var]], analysis_results$eta_results[[var]])
      }
    }
    if (!is.null(analysis_results$letters_results)) {
      if (is.null(all_results_letters[[var]])) {
        all_results_letters[[var]] <- analysis_results$letters_results[[var]]
      } else {
        all_results_letters[[var]] <- rbind(all_results_letters[[var]], analysis_results$letters_results[[var]])
      }
    }
  }
}

# Combine results for each response variable and save to separate files
for (response_var in response_vars) {
  results_kruskal_combined <- do.call(rbind, lapply(all_results_kruskal, function(x) x[which(x$measure == response_var),]))
  results_eta_combined <- do.call(rbind, lapply(all_results_eta, function(x) x[which(x$measure == response_var),]))
  results_letters_combined <- do.call(rbind, lapply(all_results_letters, function(x) x[which(x$measure == response_var),]))
  
  # Print the combined results for this response variable
  cat("Results for response variable:", response_var, "\n")
  print(results_kruskal_combined)
  print(results_eta_combined)
  print(results_letters_combined)
  
  # Write final results to files for each response variable
  write.table(results_kruskal_combined, 
              paste0(cluster_dir, "/", "kruskal_wallis_", response_var, "_by_NT_CLUSTER.txt"), 
              sep = "\t", row.names = FALSE)
  
  write.table(results_eta_combined, 
              paste0(cluster_dir, "/", "eta_squared_", response_var, "_by_NT_CLUSTER.txt"), 
              sep = "\t", row.names = FALSE)
  
  write.table(results_letters_combined, 
              paste0(cluster_dir, "/", "letters_", response_var, "_by_NT_CLUSTER.txt"), 
              sep = "\t", row.names = FALSE)
}


## FIG. S7: PLOT RATE OF CHANGE BY CLUSTER_NT (OVER THINNING LEVELS) + KRUSKAL WALLIS LETTERS (all trees vs median variables)

all_models_cluster <- read.table( file.path(cluster_dir, "RATES_CHANGES_CLUST_NT_ALL_TREES_plus_MEDIAN_MODELS.txt"), sep = "\t", header = TRUE)

all_models_cluster$clust_nt <- as.character(all_models_cluster$clust_nt)
all_models_cluster$Variable <- as.character(all_models_cluster$Variable)

letter_RCH_cluster_nt <- read.table(file.path(cluster_dir, "letters_rate_of_change_by_NT_CLUSTER.txt"), sep = "\t", header = TRUE)

letter_RCH_cluster_nt$clust_nt <- as.character(letter_RCH_cluster_nt$cluster)
letter_RCH_cluster_nt$Variable <- as.character(letter_RCH_cluster_nt$variable)

plot_data <- left_join(all_models_cluster, letter_RCH_cluster_nt, by=c("Variable","clust_nt"))

data1 <- plot_data %>% 
  group_by(Variable, clust_nt) %>% 
  summarize(RCH = max(rate_of_change, na.rm = TRUE))

# Check the column names of both dataframes
colnames(data1)
colnames(letter_RCH_cluster_nt)

# Now perform the left_join
plot_data1 <- left_join(data1, letter_RCH_cluster_nt, by = c("Variable", "clust_nt"))

# Ensure clust_nt is numeric and order the dataframe
plot_data2 <- plot_data1 %>%
  mutate(clust_nt = as.numeric(as.character(clust_nt))) %>%  # Convert clust_nt to numeric
  arrange(Variable,clust_nt)  # Arrange by clust_nt

# Check the result
head(plot_data2)

# Define the levels of Variable for which you want to perform the calculations
variable_levels <-c("F1_cbh", "F1_depth", "F1_LAD","F1_Hdepth","MXLAD_cbh", "BP_cbh", "BR_cbh",
                    "MXLAD_depth", "BP_depth", "MXLAD_dist", "BP_dist",
                    "MXLAD_Hdepth","BP_Hdepth", "MAX_height", "total_LAI", "understory_LAI")


# Convert the variable column to a factor with specified levels
plot_data2$Variable <- factor(plot_data2$Variable, levels = variable_levels)


# Convert clust_nt to a factor with specified levels
plot_data2$clust_nt <- factor(plot_data2$clust_nt, levels = c(1, 2, 3, 4))
plot_data$clust_nt <- factor(plot_data$clust_nt, levels = c(1, 2, 3, 4))


summary_data <- plot_data %>%
  group_by(Variable, clust_nt, letter) %>%  # Group by letter as well
  summarise(
    trees_change = median(rate_of_change, na.rm = TRUE),
    sd_change = sd(rate_of_change, na.rm = TRUE),
    lower_bound = quantile(rate_of_change, 0.25, na.rm = TRUE),
    upper_bound = quantile(rate_of_change, 0.75, na.rm = TRUE),
    median_change = median(rate_of_change_median, na.rm = TRUE)
  ) %>%
  mutate(
    Variable_numeric = match(Variable, responses)
  ) %>%
  ungroup()  # Ungroup to avoid grouped output warnings



# Filter out rows with NaN in lower_bound or upper_bound and remove NA Variable
summary_data <- summary_data %>%
  filter(!is.nan(lower_bound) & !is.nan(upper_bound)) %>%
  filter(complete.cases(Variable)) %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

# Create xmin and xmax for the ribbons
summary_data <- summary_data %>%
  mutate(
    xmin = Variable_numeric - 0.2,  # Adjust these values to change ribbon width
    xmax = Variable_numeric + 0.2
  )

summary_data$clust_nt <- factor(summary_data$clust_nt, levels = c(1, 2, 3, 4))
# Define cluster colors
cluster_colors <- c("1" = "orange", "2" = "brown", "3" = "green", "4" = "lightblue")

# Ensure the names of the color vector match factor levels
names(cluster_colors) <- levels(summary_data$clust_nt)


responses <-  c("F1_cbh", "F1_depth", "F1_LAD","F1_Hdepth","MXLAD_cbh", "BP_cbh", "BR_cbh",
                "MXLAD_depth", "BP_depth", "MXLAD_dist", "BP_dist",
                "MXLAD_Hdepth","BP_Hdepth", "MAX_height", "total_LAI", "understory_LAI")


# Convert Variable to a factor based on responses for correct ordering in plots
summary_data$Variable <- factor(summary_data$Variable, levels = responses)


# Fix: Use scale_x_discrete() instead of scale_x_continuous()
combined_plot <- ggplot(summary_data) +
  
  # Ribbon for upper and lower bounds for each cluster
  geom_rect(aes(xmin = as.numeric(clust_nt) - 0.4, xmax = as.numeric(clust_nt) + 0.4, 
                ymin = lower_bound, ymax = upper_bound, fill = clust_nt), 
            alpha = 0.5) +  # Increase alpha for better visibility
  
  # Points for trees_change (filled circle)
  geom_point(aes(x = clust_nt, y = trees_change, shape = "Trees Change"), 
             size = 3, color = "darkgreen") +
  
  # Points for median_change (asterisks)
  geom_point(aes(x = clust_nt, y = median_change, shape = "Median Change"), 
             size = 4, color = "darkgreen") +
  
  # Add a dashed line at y = 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.7) +
  
  # Plot labels
  labs(x = "Cluster NT",
       y = "Median Rate of Change (RCH)",
       title = "Median Rate of Change by Cluster") +
  
  # X-axis breaks and labels (Fix applied)
  scale_x_discrete() +  # Use discrete scale for factor variables
  
  # Use cluster_colors for the fill of the ribbon
  scale_fill_manual(values = cluster_colors, guide = "none") +  # Remove fill legend
  
  # Remove color legends
  guides(color = "none") +
  
  # Custom shapes for legend (only for change types)
  scale_shape_manual(values = c("Trees Change" = 16, "Median Change" = 8), name = "Change Type") +
  
  # Add significance letters slightly above the error bars with a 5% increase
  geom_text(aes(x = clust_nt, y = upper_bound - (upper_bound *0.05), label = letter), 
            vjust = -0.2, 
            size = 4, 
            color = "black",
            nudge_y = 2.5) +  # Adjust nudge_y as needed
  
  # Faceting by Variable
  facet_wrap(~ Variable, scales = "free_y", 
             labeller = labeller(Variable = as_labeller(setNames(responses, responses)))) +
  
  # Styling
  theme_minimal() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, size = 14),
    axis.text.y = element_text(angle = 0, size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 14, color = "black"),  # Increase size of facet labels
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 16, face = "bold"),  # Increase legend title size
    legend.position = "right"  # Adjust legend position as needed
  )

# Print the combined plot
print(combined_plot)


# Save the combined plot as one image
ggsave(filename = file.path(cluster_dir, "FIG_S6_COMPARISON_RATES_CHANGE_CLUST_NT_p2575_ALL_TREES_vs_MEDIAN_VALUES.png"), 
       plot = combined_plot, width = 12, height = 12, dpi=300)






