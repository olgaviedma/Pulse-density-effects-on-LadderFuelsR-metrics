
#**Impact of LiDAR Pulse Density on Forest Fuels Metrics derived using LadderFuelsR**

#**Authors: Olga Viedma and JM Moreno**

#This study evaluates the impact of LiDAR pulse density thinning on forest structure characterization using the LadderFuelsR framework. 
#High-resolution LiDAR data were collected from Mediterranean forest sites with diverse structures and fire histories and systematically thinned to simulate varying pulse densities. 
#Key forest metrics, including leaf area density (LAD), leaf area index (LAI), canopy base height, depth, height of fuel layers, and inter-layer distances, were derived to assess changes at different thinning levels.
#This sixth section deals with the estimation of the rate of change [RCH] of the main fuel properties from each individual tree. 
#Regression models were applied to each tree individually (over thinned levels, n=10) to get the median RCH from the best models.

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

################################################################
#28. FITTED MODELS FOR EACH TREE (using all thinned levels for each tree)
################################################################

# Load necessary libraries
library(dplyr)
library(broom)

metrics_crowns_dir <- "RESULTS/MODELS"

metrics_file1 <- read.table("WATER_LAI/ALL_METRICS_AND_DIFFERENCES_NO_THINNED_ALL_SITES_depurated.txt", sep = "\t", header = TRUE)

table(metrics_file1$sites, metrics_file1$thin_level)

metrics_file1$thin_level<-factor( metrics_file1$thin_level)
metrics_file1$sites<-factor( metrics_file1$sites)
metrics_file1$treeID1<-factor( metrics_file1$treeID1)
metrics_file1$treeID<-factor( metrics_file1$treeID)


metrics_file2 <- metrics_file1[, c("treeID","treeID1","thin_level" ,"sites", "Hcbh1","dptf1","Hcb1_H1","Hdptf1",  "mxld_Hc","bp_Hcbh","Hcbh_br","mxld_dp","bp_dptf","mxld_ff","bp_ffds", "mxld_Hdp", "bp_Hdpt","mx_hgh_", "lai_tot", "undrst_")]

names(metrics_file2) <- c("treeID","treeID1","thin_level" ,"sites","F1_cbh", "F1_depth", "F1_LAD","F1_Hdepth","MXLAD_cbh", "BP_cbh", "BR_cbh", "MXLAD_depth", "BP_depth", "MXLAD_dist", "BP_dist", "MXLAD_Hdepth","BP_Hdepth", "MAX_height", "total_LAI", "understory_LAI")

# Define the function to reorder levels
reorder_levels <- function(group, thin_level_order) {
  group$thin_level <- factor(group$thin_level, levels = thin_level_order)
  return(group)
}

# Define the order of thin_level levels
thin_level_order <- c("NO_THINNED",
                      "THINNED_100P",
                      "THINNED_50P",
                      "THINNED_25P",
                      "THINNED_10P",
                      "THINNED_5P",
                      "THINNED_4P",
                      "THINNED_3P",
                      "THINNED_2P",
                      "THINNED_1P")


# Initialize lists to store results
all_exp_models <- list()
all_linear_models <- list()
all_log_models <- list()

# Get unique site names
sites_list <- unique(metrics_file2$sites)

# Loop through each site
for (site in sites_list) {
  # Filter the dataset for the current site
  site_data <- subset(metrics_file2, sites == site)
  
  # Loop through each treeID in the current site
  for (treeID_level in unique(site_data$treeID1)) {
    # Filter the dataset for the current treeID
    subset_df <- site_data[site_data$treeID1 == treeID_level, ]
    
    # Reorder levels of thin_level variable
    subset_df <- subset_df[complete.cases(subset_df$thin_level), ]
    subset_df$thin_level <- factor(subset_df$thin_level, levels = thin_level_order)
    subset_df$thin_level <- droplevels(subset_df$thin_level)
    subset_df <- subset_df[order(subset_df$thin_level), ]
    
    # Ensure the treeID and site are valid
    current_treeID <- treeID_level
    current_site <- site
    
    # Drop unnecessary columns
    subset_df <- subset_df[, !names(subset_df) %in% c("treeID","treeID1", "thin_level", "sites")]
    
    # Iterate over each column name in the dataset
    for (col_name in names(subset_df)) {
      col <- subset_df[[col_name]]
      
      # Convert factor variables to numeric if applicable
      if (is.factor(col)) {
        col <- as.numeric(as.character(col))
      }
      
      # Check if the variable is numeric and has the required length
      if (is.numeric(col) && length(col) == 10 && all(!is.na(col) & is.finite(col))) {
        # Check if there are any non-zero elements
        if (any(col != 0)) {
          # Create an index vector as the independent variable
          x <- seq_along(col)
          
          ### Exponential Model
          start_a <- max(col, na.rm = TRUE)
          start_b <- 0.1
          exp_model <- try(nls(col ~ a * exp(b * x), start = list(a = start_a, b = start_b)), silent = TRUE)
          if (!inherits(exp_model, "try-error")) {
            SS_tot <- sum((col - mean(col))^2)
            SS_res <- sum((col - predict(exp_model))^2)
            R_squared <- 1 - (SS_res / SS_tot)
            coefficients <- coef(exp_model)
            all_exp_models[[length(all_exp_models) + 1]] <- data.frame(
              Variable = col_name,
              Slope = coefficients[2],
              Intercept = coefficients[1],
              R_squared = R_squared,
              treeID = current_treeID,
              sites = current_site
            )
          }
          
          ### Logarithmic Model
          start_a <- min(col, na.rm = TRUE)
          start_b <- 1
          log_model <- try(nls(col ~ a + b * log(x), start = list(a = start_a, b = start_b)), silent = TRUE)
          if (!inherits(log_model, "try-error")) {
            SS_tot <- sum((col - mean(col))^2)
            SS_res <- sum((col - predict(log_model))^2)
            R_squared <- 1 - (SS_res / SS_tot)
            coefficients <- coef(log_model)
            all_log_models[[length(all_log_models) + 1]] <- data.frame(
              Variable = col_name,
              Slope = coefficients[2],
              Intercept = coefficients[1],
              R_squared = R_squared,
              treeID = current_treeID,
              sites = current_site
            )
          }
          
          ### Linear Model
          linear_model <- lm(col ~ x)
          if (!inherits(linear_model, "try-error")) {
            SS_tot <- sum((col - mean(col))^2)
            SS_res <- sum((col - predict(linear_model))^2)
            R_squared <- 1 - (SS_res / SS_tot)
            coefficients <- coef(linear_model)
            all_linear_models[[length(all_linear_models) + 1]] <- data.frame(
              Variable = col_name,
              Slope = coefficients[2],
              Intercept = coefficients[1],
              R_squared = R_squared,
              treeID = current_treeID,
              sites = current_site
            )
          }
        }
      }
    }
  }
}


# Combine all results into data frames
results_exp <- do.call(rbind, all_exp_models)
results_exp$sites <- factor(results_exp$sites)
results_exp$Variable <- factor(results_exp$Variable)
results_exp$model <- rep("exp", nrow(results_exp))

results_linear <- do.call(rbind, all_linear_models)
results_linear$model <- rep("lin", nrow(results_linear))

results_log <- do.call(rbind, all_log_models)
results_log$model <- rep("log", nrow(results_log))

# Write summary results to file
write.table(results_linear, file.path(metrics_crowns_dir,"FITTED_LIN_FUNCTIONS_ALL_TREES.txt"), sep = "\t", row.names = FALSE)
write.table(results_exp, file.path(metrics_crowns_dir, "FITTED_EXP_FUNCTIONS_ALL_TREES.txt"), sep = "\t", row.names = FALSE)
write.table(results_log, file.path(metrics_crowns_dir, "FITTED_LOG_FUNCTIONS_ALL_TREES.txt"), sep = "\t", row.names = FALSE)

table(results_linear$Variable <-factor(results_linear$Variable))
table(results_exp$Variable <-factor(results_exp$Variable))
table(results_log$Variable <-factor(results_log$Variable))



################################################################
#29.JOINING ALL FITTED MODELS FOR EACH TREE AND EXTRACT THE BEST R2
################################################################

### JOINING all FIT MODELS AND EXTRACT THE BEST MODEL FOR EACH VARIABLE AND TREE

crowns_dir <- "RESULTS/MODELS"

# Read the models dataframes
models_list <- list.files(crowns_dir, pattern = glob2rx("FITTED_*_FUNCTIONS_ALL_TREES.txt"), full.names = TRUE)
selected_models <- grep("EXP|LIN|LOG", models_list, value = TRUE)

# Create an empty list to store processed dataframes
output_list <- list()

# Loop over each file in models_list
for (i in seq_along(selected_models)) { 
  
  # Read data from the file
  data <- read.table(selected_models[i], sep = "\t", header = TRUE)
  
  cleaned_data <- data[!is.na(data$R_squared) & data$R_squared != -Inf, ]
  
  # Subset the dataframe to include only the desired columns
  data1 <- cleaned_data %>%
    dplyr::select(treeID, sites, Variable,Slope,Intercept, R_squared, model) %>%
    # Convert non-numeric R_squared values to NA
    mutate(R_squared = as.numeric(as.character(R_squared))) %>%
    # Replace NA values with 1.0 (or any other default value)
    mutate(R_squared = ifelse(is.na(R_squared), 0.0, R_squared)) %>%
    # Cap values to be between 0.0 and 1.0
    mutate(R_squared = pmax(0.0, pmin(1.0, R_squared)))
  
  
  # Append the subset dataframe to output_list
  output_list[[i]] <- data1
}

# Combine all dataframes into one
combined_df_all <- do.call(rbind, output_list)
combined_df_all$treeID<-factor(combined_df_all$treeID)
combined_df_all$sites<-factor(combined_df_all$sites)
combined_df_all$Variable<-factor(combined_df_all$Variable)
combined_df_all$model<-factor(combined_df_all$model)

write.table(combined_df_all, file.path(crowns_dir, "R2_FITTED_ALL_FUNCTIONS_ALL_TREES.txt"), sep = "\t", row.names = FALSE)


# Find the best model for each variable

library(dplyr)

# Select the best model (max R_squared) for each tree, variable, and site
best_models <- combined_df_all %>%
  group_by(Variable, treeID, sites) %>%
  slice_max(R_squared, with_ties = FALSE) %>%  # Retain only the single best result
  arrange(Variable, treeID, sites)

# Write the output to a file
write.table(best_models, file.path(crowns_dir, "BEST_R2_FITTED_ALL_FUNCTIONS_ALL_TREES.txt"), sep = "\t", row.names = FALSE)


################################################################
#30.RATES OF CHANGE of EACH VARIABLE FOR EACH TREE (SLOPE from the BEST MODEL) (TABLE_S5 & FIG S5)
################################################################

library(dplyr)

crowns_dir <- "RESULTS/MODELS"

best_models<- read.table (file.path(crowns_dir, "BEST_R2_FITTED_ALL_FUNCTIONS_ALL_TREES.txt"), sep = "\t", header = TRUE)
best_models$treeID<-factor(best_models$treeID)
best_models$sites<-factor(best_models$sites)
best_models$Variable<-factor(best_models$Variable)
best_models$model<-factor(best_models$model)


# Define the desired order of responses
responses <-c("F1_cbh", "F1_depth", "F1_LAD","F1_Hdepth","MXLAD_cbh", "BP_cbh", "BR_cbh",
              "MXLAD_depth", "BP_depth", "MXLAD_dist", "BP_dist",
              "MXLAD_Hdepth","BP_Hdepth", "MAX_height", "total_LAI", "understory_LAI")


best_models1 <- best_models %>%
  dplyr::filter(Variable %in% responses)

# Reorder factor levels for the Variable column based on the desired order
best_models1$Variable <- factor(best_models1$Variable, levels = responses)

# Reorder the rows based on 'Variable'
best_models1 <- best_models1[order(best_models1$Variable), ]



# RATES OF CHANGE FROM ALL TREES

# Function to calculate percentage change for exponential models
calculate_percentage_change_exp <- function(Slope) {
  return((exp(Slope) - 1) * 100)
}

# Function to calculate the rate of change for logarithmic models
calculate_rate_of_change_log <- function(Slope) {
  return(Slope * 0.01)  # Slope gives the percentage change in response to a 1% change in the predictor
}

# Updated best_models_clean1 calculation with rate of change for different models
best_models_change <- best_models1 %>%
  mutate(rate_of_change = case_when(
    model == "lin" ~ Slope * 100,  # Linear model rate of change
    model == "exp" ~ calculate_percentage_change_exp(Slope),  # Exponential model rate of change
    model == "log" ~ calculate_rate_of_change_log(Slope),  # Logarithmic model rate of change
    TRUE ~ NA_real_  # Default to NA for other models
  )) %>%
  ungroup()  # Ungroup after the calculation

# Print the updated table with percentage change
print(best_models_change)

write.table(best_models_change, file.path(crowns_dir, "RATES_CHANGES_BEST_R2_ALL_TREES.txt"), sep = "\t", row.names = FALSE)



# TABLE S5. CONVERT BEST MODELS for ALL TREES AND RATES OF CHANGE (SLOPE) INTO A .DOCX TABLE

library(officer)
library(flextable)

crowns_dir <- "RESULTS/MODELS"

reference_values_wide<- read.table (file.path(crowns_dir, "MEDIAN_VARIABLES_BY_TREES_quantiles0595_WIDE_all_sites.txt"), sep = "\t", header = TRUE)

reference_values_wide$Variable<-factor(reference_values_wide$Variable)
reference_values_wide$treeID<-factor(reference_values_wide$treeID)
reference_values_wide$sites<-factor(reference_values_wide$sites)

best_models_change<- read.table (file.path(crowns_dir, "RATES_CHANGES_BEST_R2_ALL_TREES.txt"), sep = "\t", header = TRUE)

best_models_change$Variable<-factor(best_models_change$Variable)
best_models_change$treeID<-factor(best_models_change$treeID)
best_models_change$sites<-factor(best_models_change$sites)

best_models_change2 <- best_models_change %>%
  left_join(reference_values_wide, by = c( "treeID","Variable", "sites"))


# Create a Word document
doc1 <- read_docx()

# Define a helper function for get_mode if not already defined
get_mode <- function(x) {
  uniq_x <- unique(x)
  uniq_x[which.max(tabulate(match(x, uniq_x)))]
}

# Updated pipeline with corrections
summary_data3 <- best_models_change2 %>%
  dplyr::group_by(Variable) %>%
  dplyr::summarize(
    rate_change_trees1 = median(rate_of_change, na.rm = TRUE),
    Slope_trees1 = median(Slope, na.rm = TRUE),
    sd_slope_trees1 = sd(Slope, na.rm = TRUE),
    Intercept_trees1 = median(Intercept, na.rm = TRUE),
    sd_intercept_trees1 = sd(Intercept, na.rm = TRUE),
    R2_trees1 = median(R_squared, na.rm = TRUE),
    sd_R2_trees1 = sd(R_squared, na.rm = TRUE),
    Q1_R2 = quantile(R_squared, 0.25, na.rm = TRUE),
    Q3_R2 = quantile(R_squared, 0.75, na.rm = TRUE),
    mediana1 = median(median, na.rm = TRUE),
    sd_mediana1 = median(sd, na.rm = TRUE),
    model1 = get_mode(model),  # Use the most frequent model for each Variable
    # Calculate Q1 and Q3 for rate_of_change
    Q1_Rtch = quantile(rate_of_change, 0.25, na.rm = TRUE),
    Q3_Rtch = quantile(rate_of_change, 0.75, na.rm = TRUE)
  ) %>%
  mutate(
    # Calculate bounds using Q1 and Q3
    lower_bound = Q1_Rtch,
    upper_bound = Q3_Rtch,
    lower_bound_R2 = Q1_R2,
    upper_bound_R2 = Q3_R2,
    # Convert Variable to a numeric factor
    Variable_numeric = as.numeric(factor(Variable))
  )

# Check the output
print(summary_data3)
summary_data3 <- summary_data3 %>% 
  mutate(across(where(is.numeric), ~ formatC(., format = "f", digits = 2)))

summary_data4 <- summary_data3 %>%
  dplyr::mutate( Intercept = paste0(Intercept_trees1, " ± ", sd_intercept_trees1),
                 Slope = paste0(Slope_trees1, " ± ", sd_slope_trees1),
                 R2 = R2_trees1,
                 bounds_R2 = paste0("(", lower_bound_R2, ",", upper_bound_R2, ")"),
                 #rate_change_trees = paste0(rate_change_trees1, " ± ", sd_change_trees1),
                 RCH = paste0(rate_change_trees1),
                 bounds_rate = paste0("(", lower_bound, ",", upper_bound, ")"),
                 median = paste0(mediana1, " ± ", sd_mediana1),
                 model= model1) %>%  # closed parentheses correctly
  arrange(Variable) %>%
  dplyr::select(Variable,  model, Intercept, Slope, R2, bounds_R2, RCH,bounds_rate, median)


summary_data4$Variable<- factor(summary_data4$Variable, levels=responses)

summary_data4 <- summary_data4[order(summary_data4$Variable), ]

# Convert the tibble to a flextable
ft1 <- flextable(summary_data4)
autofit(ft1)  # Automatically adjusts the column width to fit the content

# Add the flextable to the document
doc1 <- doc1 %>% 
  body_add_flextable(ft1)

# Save the document
print(doc1, target = file.path(crowns_dir,"TABLE_S5_table_best_models_RATE_CHANGE_P2575_using_ALL_TREES.docx"))



# FIG. S5: PLOTS OF COMPARISON BETWEEN RATES OF CHANGE (BAR PLOTS) FROM MODELS OF EACH ALL TREES AND FROM FITTED MODELS OVER MEDIANS DIRECTLY

library(ggplot2)
library(RColorBrewer)
library(viridis)  # Load the viridis package for color generation
library(scales)
library(ggsci)

best_models_medians<- read.table (file.path(crowns_dir, "RATES_CHANGES_BEST_R2_MEDIAN_VALUES_all_sites.txt"), sep = "\t", header = TRUE)
best_models_medians$Variable<-factor(best_models_medians$Variable)
best_models_medians$model<-factor(best_models_medians$model)
best_models_medians <- best_models_medians %>%rename_with(~ paste0(., "_median"),  .cols = -c(Variable, model))

best_models_change <- read.table(file.path(crowns_dir,"RATES_CHANGES_BEST_R2_ALL_TREES.txt"), sep = "\t", header = TRUE)
best_models_change$Variable<-factor(best_models_change$Variable)
best_models_change$model<-factor(best_models_change$model)

# Define the desired order of responses
responses <- c("F1_cbh", "F1_depth", "F1_LAD","F1_Hdepth","MXLAD_cbh", "BP_cbh", "BR_cbh",
               "MXLAD_depth", "BP_depth", "MXLAD_dist", "BP_dist",
               "MXLAD_Hdepth", "BP_Hdepth","MAX_height", "total_LAI", "understory_LAI")

# Set factor levels
best_models_change$Variable <- factor(best_models_change$Variable, levels = responses)
best_models_medians$Variable <- factor(best_models_medians$Variable, levels = responses)

best_models_change <- na.omit(best_models_change)
best_models_medians <- na.omit(best_models_medians)

# Check the summary to confirm changes
summary(best_models_change)

# Order the dataframe by the reordered Variable factor
best_models_change <- best_models_change[order(best_models_change$Variable), ]
best_models_medians <- best_models_medians[order(best_models_medians$Variable), ]

best_models_trees_median <- dplyr::left_join(best_models_change, best_models_medians, by = c( "Variable"))

# Calculate summary statistics, including quantiles
summary_data <- best_models_trees_median  %>%
  group_by(Variable) %>%
  summarise(
    trees_change = median(rate_of_change, na.rm = TRUE),
    sd_change = sd(rate_of_change, na.rm = TRUE),
    lower_bound = quantile(rate_of_change, 0.25, na.rm = TRUE),  # Corrected here
    upper_bound = quantile(rate_of_change, 0.75, na.rm = TRUE),  # Corrected here
    median= median(median_median, na.rm = TRUE),
    median_change = median(rate_of_change_median, na.rm = TRUE)
  ) %>%
  mutate(
    Variable_numeric = match(Variable, responses)  # Assign numeric order based on 'responses' vector
  )


# Filter out rows with NaN in lower_bound or upper_bound and remove NA Variable
summary_data <- summary_data %>%
  filter(!is.nan(lower_bound) & !is.nan(upper_bound)) %>%
  filter(complete.cases(Variable)) %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

write.table(summary_data, file.path(crowns_dir, "MEDIAN_RATES_CHANGES_quantile2575_ALL_TREES_MEDIAN.txt"), sep = "\t", row.names = FALSE)

# Create xmin and xmax for the ribbons
summary_data <- summary_data %>%
  mutate(
    xmin = Variable_numeric - 0.2,  # Adjust these values to change ribbon width
    xmax = Variable_numeric + 0.2
  )


# Ensure Variable is treated as character to avoid issues with factors
summary_data$Variable <- as.character(summary_data$Variable)


responses <- c("F1_cbh", "F1_depth", "F1_LAD","F1_Hdepth","MXLAD_cbh", "BP_cbh", "BR_cbh",
               "MXLAD_depth", "BP_depth", "MXLAD_dist", "BP_dist",
               "MXLAD_Hdepth", "BP_Hdepth","MAX_height", "total_LAI", "understory_LAI")


summary_data$Variable <- factor(summary_data$Variable, levels = responses, ordered = TRUE)

summary_data <- summary_data %>%
  arrange(Variable) %>%
  mutate(Variable_numeric = as.numeric(Variable))


# Get unique variable names from your summary data
unique_variables <- unique(summary_data$Variable)

num_unique <- length(unique_variables)


variable_colors <- scales::hue_pal()(length(unique_variables))

# Assign colors to variables dynamically
names(variable_colors) <- unique_variables

p <- ggplot() +
  # Ribbon using geom_rect() - fill should not add to legend for points
  geom_rect(data = summary_data, aes(xmin = Variable_numeric - 0.4, xmax = Variable_numeric + 0.4, 
                                     ymin = lower_bound, ymax = upper_bound, fill = Variable), 
            alpha = 0.5) +  # Increase alpha for better visibility
  
  # Points for median changes (filled circle) with shape mapped for legend
  geom_point(data = summary_data, aes(x = Variable_numeric, y = trees_change, shape = "Trees Change"), 
             size = 3, color = "darkgreen") +
  
  # Points for median_change as asterisks
  geom_point(data = summary_data, aes(x = Variable_numeric, y = median_change, shape = "Median Change"), 
             size = 4, color = "darkgreen") +  # Shape 8 for asterisk
  
  labs(title = "Percentage Change by Variable with Bounds",
       x = "Variable",
       y = "Median Percentage Change") +
  
  # Set x breaks and labels directly from the summary data, ensuring alignment
  scale_x_continuous(breaks = summary_data$Variable_numeric, 
                     labels = summary_data$Variable, 
                     expand = c(0, 0)) +  # Remove expansion to align labels exactly with points
  
  # Custom fill colors for the geom_rect()
  scale_fill_manual(values = variable_colors, name = "Variable") +
  
  # Custom shapes for legend
  scale_shape_manual(values = c("Trees Change" = 16, "Median Change" = 8), name = "Change Type") +
  
  # Remove fill scale for the geom_rect() to avoid variable names in legend
  guides(fill = "none") +  # This line ensures the fill legend is suppressed
  
  theme_minimal() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 14),  # Center align labels at marks
    axis.text.y = element_text(angle = 0, size = 14),
    axis.title.x = element_text(size = 14),                        # Adjust size of x-axis title
    axis.title.y = element_text(size = 14),                        # Adjust size of y-axis title
    plot.title = element_text(hjust = 0.5),
    legend.position = "right" ,  # Increase legend title size
    legend.text = element_text(size = 14)    # Increase legend text size
  )

# Print the plot
print(p)

# Save the plot
ggsave(filename = file.path(crowns_dir, paste0("FIG_S5_PLOT_COMPARISON_RATES_CHANGE_TREES_MEDIAN_VAR_p2575_THIN_LOWEST.png")), plot = p, width = 8, height = 6)



## PLOT RATES OF CHANGE (MEDIAN CURVES) BASED ON FITTED MODELS OVER INDIVIDUAL TREES

library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(ggrepel)

crowns_dir <- "RESULTS/MODELS"

all_data_lidarpod1 <- read.table(file.path(crowns_dir,"ALL_METRICS_together_long_NO_THINNED_depurated_names.txt"), sep = "\t", header = TRUE)


all_data_lidarpod1$treeID <- factor(all_data_lidarpod1$treeID,
                                    levels=unique(all_data_lidarpod1$treeID))
all_data_lidarpod1$thin_level <- factor(all_data_lidarpod1$thin_level,
                                        levels=unique(all_data_lidarpod1$thin_level))
all_data_lidarpod1$sites <- factor(all_data_lidarpod1$sites,
                                   levels=unique(all_data_lidarpod1$sites))

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
all_data_lidarpod1$thin_level <- factor(all_data_lidarpod1$thin_level, levels = new_order)

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
levels(all_data_lidarpod1$thin_level) <- new_names

responses_def <- c("F1_cbh", "F1_depth", "F1_LAD","F1_Hdepth","MXLAD_cbh", "BP_cbh", "BR_cbh",
                   "MXLAD_depth", "BP_depth", "MXLAD_dist", "BP_dist",
                   "MXLAD_Hdepth", "BP_Hdepth", "MAX_height",
                   "total_LAI", "understory_LAI")

numeric_data <- all_data_lidarpod1 %>%
  dplyr::select(thin_level, where(is.numeric)) %>%
  dplyr::select(-treeID1)  # Exclude treeID1 explicitly


########### RATES OF CHANGE AND R2 VALUES ##################3

best_models_trees2 <- read.table(file.path(crowns_dir, "RATES_CHANGES_BEST_R2_ALL_TREES.txt"), sep = "\t", header = TRUE)


numeric_data_models <- best_models_trees2 %>%
  dplyr::select(Variable, where(is.numeric)) %>%
  dplyr::select(-treeID)  # Exclude treeID1 explicitly

# Create filtered data for the current variable
filtered_data_long_models <- numeric_data_models %>%
  group_by(Variable) %>%
  summarize(across(everything(), 
                   list(median = ~ median(.x, na.rm = TRUE),
                        sd = ~ sd(.x, na.rm = TRUE),
                        lower_ci = ~ quantile(.x, 0.05, na.rm = TRUE),
                        upper_ci = ~ quantile(.x, 0.95, na.rm = TRUE)),
                   .names = "{.col}_{.fn}"),
            .groups = 'drop') 

#Define the desired order of responses
responses <-c("F1_cbh", "F1_depth", "F1_LAD","F1_Hdepth","MXLAD_cbh", "BP_cbh", "BR_cbh",
              "MXLAD_depth", "BP_depth", "MXLAD_dist", "BP_dist",
              "MXLAD_Hdepth","BP_Hdepth", "MAX_height", "total_LAI", "understory_LAI")

# Reorder factor levels for the Variable column based on the desired order
filtered_data_long_models$Variable <- factor(filtered_data_long_models$Variable, levels = responses)

# Reorder the rows based on 'Variable'
filtered_data_long_models <- filtered_data_long_models[order(filtered_data_long_models$Variable), ]



########  PLOTTING  ###############3

# Initialize a list to hold all individual plots across all groups
all_plots <- list()

# Set the output directory for saving plots
output_dir <- "RESULTS/MODELS"

# Initialize a dataframe to hold data for the faceted plot
all_variable_data <- data.frame()
annotations_list <- data.frame()  # Initialize a dataframe to hold annotations

# Loop over each variable in responses
for (var in responses_def) {
  
  # Create filtered data for the current variable
  filtered_data_long <- numeric_data %>%
    dplyr::select(all_of(c(var, "thin_level"))) %>%
    group_by(thin_level) %>%
    summarize(across(everything(), 
                     list(media = ~ median(.x, na.rm = TRUE),
                          sd = ~ sd(.x, na.rm = TRUE),
                          lower_ci = ~ quantile(.x, 0.05, na.rm = TRUE),
                          upper_ci = ~ quantile(.x, 0.95, na.rm = TRUE)),
                     .names = "{.col}_{.fn}"),
              .groups = 'drop') %>%
    pivot_longer(
      cols = -thin_level,
      names_to = c("variable", "stat"),
      names_pattern = "(.*)_(media|sd|lower_ci|upper_ci)"
    ) %>%
    pivot_wider(
      names_from = stat,
      values_from = value,
      names_glue = "{stat}"
    ) %>%
    mutate(
      lower_bound = lower_ci,
      upper_bound = upper_ci
    )
  
  # Set the order of the variable factor to match the responses vector
  #filtered_data_long$variable <- factor(filtered_data_long$variable, levels = responses)
  
  
  # Append the current filtered data to the combined dataframe
  all_variable_data <- rbind(all_variable_data, filtered_data_long)
  
  # Check if there is data to plot
  if (nrow(filtered_data_long) == 0) {
    message(paste("No data for variable:", var))
    next  # Skip to the next variable if there's no data
  }
  
  # Create the plot
  plot <- ggplot(filtered_data_long, aes(x = thin_level, y = media, color = variable, group = 1)) +
    geom_line(size = 1.5, alpha = 1) +  
    geom_point(size = 3, alpha = 1) +   
    geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound, fill = variable), alpha = 0.1) +
    labs(title = paste("Variable:", var), 
         x = "thin_level", y = "Value") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )
  
  # Extract R-squared and rate of change values for the current variable
  rsq_values <- filtered_data_long_models %>%
    dplyr::filter(Variable == var) %>%
    dplyr::select(Variable, rate_of_change_median, R_squared_median) %>%
    distinct() %>%
    mutate(rate_of_change = as.numeric(rate_of_change_median),
           R_squared = as.numeric(R_squared_median))
  
  # Check if R-squared values are available
  if (nrow(rsq_values) > 0) {
    # Calculate the first thin_level and value for annotation
    first_value <- filtered_data_long %>%
      summarize(first_value = first(media[!is.na(media)]), 
                first_thin_level = first(thin_level), .groups = 'drop')
    
    # Create a new data frame for annotation
    new_annotation <- data.frame(
      thin_level = factor(first_value$first_thin_level), 
      value_at_first_thin_level = first_value$first_value, 
      R_squared = round(rsq_values$R_squared[1], 2),  # Use the first R² value
      rate_of_change = round(rsq_values$rate_of_change[1], 2),  # Use the first rate of change
      variable = var
    )
    
    # Append annotations to the annotations_list
    annotations_list <- rbind(annotations_list, new_annotation)
    
    # Add annotations to the plot
    plot <- plot +
      geom_text_repel(data = new_annotation, 
                      aes(x = thin_level, 
                          y = value_at_first_thin_level, 
                          label = paste0("R²: ", R_squared, "\nRate Change: ", rate_of_change)),  # Use new line for clarity
                      size = 3, color = "blue",
                      nudge_y = 2, 
                      box.padding = 2,
                      point.padding = 2,
                      segment.color = 'grey50',
                      force = 2.5, 
                      max.overlaps = 10)
  }
  
  # Save each plot to a unique file
  #ggsave(filename = file.path(output_dir, paste0("plot_", var, ".png")), plot = plot, width = 8, height = 6)
  
  # Store the individual plot in the list
  all_plots[[var]] <- plot
}

# Create a faceted plot for all variables

# Set factor levels in both datasets explicitly to match `responses1`
all_variable_data$variable <- factor(all_variable_data$variable, levels = responses_def)
annotations_list$variable <- factor(annotations_list$variable, levels = responses_def)


# Faceted Plot with Annotations
faceted_plot <- ggplot(all_variable_data, aes(x = thin_level, y = media, color = variable, group = variable)) +
  geom_line(size = 1.5, alpha = 0.7) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound, fill = variable, group = variable), alpha = 0.3) +
  labs(title = "", x = "Thinning levels", y = "Rate of change % (RCH)") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    legend.position = "none",
    plot.title = element_text(size = 14, hjust = 0.5, color = "black"),
    strip.text = element_text(size = 14, color = "black"),
    plot.margin = margin(10, 10, 10, 10),
    panel.spacing = unit(1, "lines")
  ) +
  facet_wrap(~ variable, scales = "free_y", labeller = labeller(variable = as_labeller(setNames(responses_def, responses_def)))) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  coord_cartesian(clip = "off")  # Prevents clipping

# Add Annotations Separately
faceted_plot <- faceted_plot +
  geom_text_repel(
    data = annotations_list,
    aes(
      x = thin_level,
      y = value_at_first_thin_level,
      label = paste0("R²: ", R_squared, " RCH: ", rate_of_change)
    ),
    size = 5, color = "black",
    nudge_y = 2, 
    box.padding = 2,
    point.padding = 2,
    segment.color = 'grey50',
    force = 2.5, 
    max.overlaps = 10
  )


# Save the Faceted Plot
ggsave(
  filename = file.path(output_dir, "PLOT_RATES_CHANGE_BEST_FITTED_MODELS_ALL_TREES_p0595_MEDIAN_VALUES.png"),
  plot = faceted_plot,
  width = 12,
  height = 8,
  dpi = 300
)

