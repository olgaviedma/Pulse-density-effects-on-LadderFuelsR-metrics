
#**Impact of LiDAR Pulse Density on Forest Fuels Metrics derived using LadderFuelsR**

#**Authors: Olga Viedma and JM Moreno**

#This study evaluates the impact of LiDAR pulse density thinning on forest structure characterization using the LadderFuelsR framework. 
#High-resolution LiDAR data were collected from Mediterranean forest sites with diverse structures and fire histories and systematically thinned to simulate varying pulse densities. 
#Key forest metrics, including leaf area density (LAD), leaf area index (LAI), canopy base height, depth, height of fuel layers, and inter-layer distances, were derived to assess changes at different thinning levels.
#This fifth section deals with the estimation of the rate of change [RCH] of the main fuel properties over the thinned levels performing various regression models (exponential, logarithmic and linear) 
#using the median values at each level (n = 10 thinning levels) directly, and converted the regression slopes into percentages. 

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


###################################################
#25. FITTED MODELS ON MEDIAN VARIABLES (only by thin_levels)
###################################################

library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra) # for grid.arrange
library(changepoint)

crowns_dir <- "RESULTS/MODELS"

filtered_data1 <- read.table(file.path(crowns_dir, "MEDIAN_VARIABLES_BY_THIN_LEVELS_quantiles0595_WIDE_all_sites.txt"), sep = "\t", header = TRUE)



## FITTED MODELS 


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

# Function to calculate summary statistics
calculate_summary <- function(results) {
  if (nrow(results) == 0) {
    return(data.frame(
      greater_than_55 = NA, greater_than_75 = NA, less_than_55 = NA, less_than_75 = NA,
      percentage_greater_than_55 = NA, percentage_greater_than_75 = NA,
      percentage_less_than_55 = NA, percentage_less_than_75 = NA
    ))
  }
  greater_than_55 <- sum(results$R_squared >= 0.55)
  greater_than_75 <- sum(results$R_squared >= 0.75)
  less_than_75 <- sum(results$R_squared < 0.75)
  less_than_55 <- sum(results$R_squared < 0.55)
  total_rows <- nrow(results)
  percentage_greater_than_75 <- (greater_than_75 / total_rows) * 100
  percentage_less_than_75 <- (less_than_75 / total_rows) * 100
  percentage_greater_than_55 <- (greater_than_55 / total_rows) * 100
  percentage_less_than_55 <- (less_than_55 / total_rows) * 100
  data.frame(
    greater_than_55, greater_than_75, less_than_55, less_than_75,
    percentage_greater_than_55, percentage_greater_than_75,
    percentage_less_than_55, percentage_less_than_75
  )
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
    geom_point(aes(color = "Median Values"), size = 2) +
    geom_line(aes(y = fitted_values, color = "Fitted Values")) +
    labs(title = paste("Regression Plot for", variable, "using", model_type, "model"),
         x = "Index",
         y = "Values",
         subtitle = paste("Intercept:", round(intercept, 2), "Slope:", round(slope, 2))) +
    scale_color_manual(values = c("Median Values" = "blue", "Fitted Values" = "red")) +
    theme_minimal()
  
  return(plot)
}


# Initialize lists to store results across all sites
all_exp_results <- list()
all_log_results <- list()
all_linear_results <- list()

# Reorder levels of thin_level variable
filtered_data1 <- filtered_data1[complete.cases(filtered_data1$thin_level), ]
filtered_data1$thin_level <- factor(filtered_data1$thin_level, levels = new_order)
filtered_data1$thin_level <- droplevels(filtered_data1$thin_level)
filtered_data1 <- filtered_data1[order(filtered_data1$thin_level), ]

# Initialize lists to store results and plots for each cluster
site_exp_results <- list()
site_log_results <- list()
site_linear_results <- list()

exp_plots <- list()
log_plots <- list()
linear_plots <- list()

# Iterate over each column in the data for modeling
for (col_name in names(filtered_data1)) {
  col <- filtered_data1[[col_name]]
  
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
        site_exp_results[[length(site_exp_results) + 1]] <- data.frame(
          Variable = col_name, Slope = coefficients[2], Intercept = coefficients[1], 
          R_squared = R_squared$adjusted_r_squared,
          model = "exp"
        )
      }
      
      # Fitting logarithmic model
      log_model <- tryCatch({
        nls(col ~ a + b * log(x), start = list(a = min(col, na.rm = TRUE), b = 1))
      }, error = function(e) NULL)
      
      if (!is.null(log_model)) {
        coefficients <- coef(log_model)
        R_squared <- calculate_adjusted_r_squared(col, predict(log_model), n_predictors = 1)
        site_log_results[[length(site_log_results) + 1]] <- data.frame(
          Variable = col_name, Slope = coefficients[2], Intercept = coefficients[1], 
          R_squared = R_squared$adjusted_r_squared,
          model = "log"
        )
      }
      
      # Fitting linear model
      linear_model <- tryCatch({
        lm(col ~ x)
      }, error = function(e) NULL)
      
      if (!is.null(linear_model)) {
        coefficients <- coef(linear_model)
        R_squared <- calculate_adjusted_r_squared(col, predict(linear_model), n_predictors = 1)
        site_linear_results[[length(site_linear_results) + 1]] <- data.frame(
          Variable = col_name, Slope = coefficients[2], Intercept = coefficients[1], 
          R_squared = R_squared$adjusted_r_squared,
          model = "lin"
        )
      }
    }
  }
}

# Append model results for each variable
all_exp_results <- do.call(rbind, site_exp_results)
all_linear_results <- do.call(rbind, site_linear_results)
all_log_results <- do.call(rbind, site_log_results)


# Plot after results collection for each variable in this cluster
for (model_type in c("exp", "log", "linear")) {
  model_results <- get(paste0("site_", model_type, "_results"))
  
  if (length(model_results) > 0) {
    # Iterate over the variables in model_results for plotting
    for (result in model_results) {
      var <- result$Variable
      plot <- fit_and_plot_regression(filtered_data1, var, model_type)
      
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


# Save the plots for each model type in chunks

for (model_type in c("exp", "log", "linear")) {
  plots <- get(paste0(model_type, "_plots"))
  
  if (length(plots) > 0) {
    plot_chunks <- split(plots, ceiling(seq_along(plots) / 12))
    
    # Ensure unique filenames by including site_level and cluster_level
    for (i in seq_along(plot_chunks)) {
      facet_plot <- do.call(grid.arrange, c(plot_chunks[[i]], ncol = 3))
      
      # Save the plot with filename reflecting both site_level and cluster_level
      ggsave(
        filename = file.path(crowns_dir, 
                             paste0("facet_plot_", 
                                    "model_", model_type, "_chunk", i, ".png")),
        plot = facet_plot,
        width = 20, 
        height = 15
      )
    }
  }
}


# Combine results across all sites
combined_exp_results <- do.call(rbind, all_exp_results)
combined_log_results <- do.call(rbind, all_log_results)
combined_linear_results <- do.call(rbind, all_linear_results)

combined_exp_results_df <- as.data.frame(t(combined_exp_results))
combined_log_results_df <- as.data.frame(t(combined_log_results))
combined_lin_results_df <- as.data.frame(t(combined_linear_results))

# Write the combined results to a single file for each model
write.table(combined_lin_results_df, file.path(crowns_dir, "FITTED_LIN_FUNCTIONS_MEDIAN_VARIABLES_all_sites.txt"), sep = "\t", row.names = FALSE)
write.table(combined_exp_results_df, file.path(crowns_dir, "FITTED_EXP_FUNCTIONS_MEDIAN_VARIABLES_all_sites.txt"), sep = "\t", row.names = FALSE)
write.table(combined_log_results_df, file.path(crowns_dir, "FITTED_LOG_FUNCTIONS_MEDIAN_VARIABLES_all_sites.txt"), sep = "\t", row.names = FALSE)


###################################################
#26. JOINING ALL FITTED MODELS AND EXTRACT THE BEST R2
###################################################

library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra) # for grid.arrange
library(changepoint)


## JOINING ALL FITTED MODELS

crowns_dir <- "RESULTS/MODELS"

# Read the models dataframes
models_list <- list.files(crowns_dir, pattern = glob2rx("FITTED_*_FUNCTIONS_MEDIAN_VARIABLES_all_sites.txt"), full.names = TRUE)
models_file <- lapply(models_list, function(X) read.table(X, sep = "\t", header = TRUE))

all_models <- do.call (rbind, models_file)

write.table(all_models, file.path(crowns_dir, "ALL_FITTED_FUNCTIONS_MEDIAN_VARIABLES_all_sites.txt"), sep = "\t", row.names = FALSE)

all_models$model <-factor(all_models$model)
all_models$Variable <-factor(all_models$Variable)


#  CALCULATE BEST MODELS

# Select the best model (max R_squared) for each tree, variable, and site
best_models_median <- all_models %>%
  group_by(Variable) %>%
  slice_max(R_squared, with_ties = FALSE) %>%  # Retain only the single best result
  arrange(Variable)

#best_models_clean <- best_models %>%dplyr::filter(R_squared != -Inf)

best_models_median$Variable <- sub("_median$", "", best_models_median$Variable)
best_models_median$Variable <-factor(best_models_median$Variable )


###  JOIN MEDIAN VALUES WITH BEST MODELS RESULTS

filtered_data_sensor_wide <- read.table(file.path(crowns_dir, "MEDIAN_VARIABLES_quantiles0595_WIDE_all_sites.txt"), sep = "\t", header = TRUE)

# Perform the left join
best_models_median1 <- best_models_median %>%
  left_join(filtered_data_sensor_wide, by = "Variable")

best_models_median1$Variable <- factor(best_models_median1$Variable)
best_models_median1$model <- factor(best_models_median1$model)

#best_models_median1 <- best_models_median1[complete.cases(best_models_median1$median), ]

write.table(best_models_median1, file.path(crowns_dir, "BEST_R2_ALL_FITTED_FUNCTIONS_MEDIANA_VARIABLES_all_sites.txt"), sep = "\t", row.names = FALSE)


###################################################
#27.RATES OF CHANGE (slope of the best model) of EACH VARIABLE by THINNED LEVEL (median values)  (FIG. 8 and TABLE S4)
##################################################

library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra) # for grid.arrange
library(changepoint)

crowns_dir <- "RESULTS/MODELS"


#  CALCULATE RATES OF CHANGE (SLOPE) OF BEST MODELS

crowns_dir<- "RESULTS/MODELS"

best_models_median1 <- read.table(file.path(crowns_dir,"BEST_R2_ALL_FITTED_FUNCTIONS_MEDIANA_VARIABLES_all_sites.txt"), sep = "\t", header = TRUE)


# Define the desired order of responses
responses <-c("F1_cbh", "F1_depth", "F1_LAD","F1_Hdepth","MXLAD_cbh", "BP_cbh", "BR_cbh",
              "MXLAD_depth", "BP_depth", "MXLAD_dist", "BP_dist",
              "MXLAD_Hdepth","BP_Hdepth", "MAX_height", "total_LAI", "understory_LAI")

best_models_median1 <- best_models_median1 %>%
  dplyr::filter(Variable %in% responses)

# Reorder factor levels for the Variable column based on the desired order
best_models_median1$Variable <- factor(best_models_median1$Variable, levels = responses)

# Reorder the rows based on 'Variable'
best_models_median1 <- best_models_median1[order(best_models_median1$Variable), ]

# Function to calculate percentage change for exponential models
calculate_percentage_change_exp <- function(Slope) {
  return((exp(Slope) - 1) * 100)
}

# Function to calculate the rate of change for logarithmic models
calculate_rate_of_change_log <- function(Slope) {
  return(Slope * 0.01)  # Slope gives the percentage change in response to a 1% change in the predictor
}

# Updated best_models_clean1 calculation with rate of change for different models
best_models_median2 <- best_models_median1 %>%
  mutate(rate_of_change = case_when(
    model == "lin" ~ Slope * 100,  # Linear model rate of change
    model == "exp" ~ calculate_percentage_change_exp(Slope),  # Exponential model rate of change
    model == "log" ~ calculate_rate_of_change_log(Slope),  # Logarithmic model rate of change
    TRUE ~ NA_real_  # Default to NA for other models
  )) %>%
  ungroup()  # Ungroup after the calculation


# Round numeric variables to 2 decimal places
best_models_median2 <- best_models_median2 %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

# Convert Variable back to factor
best_models_median2$Variable <- factor(best_models_median2$Variable)

best_models_median2$Variable<- factor(best_models_median2$Variable, levels=responses)
best_models_median2 <- best_models_median2[order(best_models_median2$Variable), ]


write.table(best_models_median2, file.path(crowns_dir, "RATES_CHANGES_BEST_R2_MEDIAN_VALUES_all_sites.txt"), sep = "\t", row.names = FALSE)


#  TABLE S3: CONVERT BEST MODELS AND RATES OF CHANGE (SLOPE) INTO A .DOCX TABLE


library(officer)
library(flextable)

best_models_median2 <- read.table(file.path(crowns_dir, "RATES_CHANGES_BEST_R2_MEDIAN_VALUES_all_sites.txt"), sep = "\t", header = TRUE)

#Define the desired order of responses
responses <-c("F1_cbh", "F1_depth", "F1_LAD","F1_Hdepth","MXLAD_cbh", "BP_cbh", "BR_cbh",
              "MXLAD_depth", "BP_depth", "MXLAD_dist", "BP_dist",
              "MXLAD_Hdepth","BP_Hdepth", "MAX_height", "total_LAI", "understory_LAI")

# Reorder factor levels for the Variable column based on the desired order
best_models_median2$Variable <- factor(best_models_median2$Variable, levels = responses)

# Reorder the rows based on 'Variable'
best_models_median2 <- best_models_median2[order(best_models_median2$Variable), ]

# Create a Word document
doc <- read_docx()

best_models_median3 <- best_models_median2 %>%
  # Format all numeric columns to have 2 decimal places
  mutate(across(where(is.numeric), ~ formatC(., format = "f", digits = 2))) %>%
  # Combine 'median' and 'sd' into the 'mediana' column with formatted values
  mutate(mediana = paste0(median, " ± ", sd)) %>%
  arrange(Variable) %>%
  dplyr::select(Variable,  model, Intercept, Slope, R_squared, rate_of_change, mediana)

# Create the flextable
ft0 <- flextable(best_models_median3)

# Try autofit after confirming the table is generated correctly
ft0 <- autofit(ft0)

# Add the flextable to the document
doc <- doc %>% 
  body_add_flextable(ft0)

# Save the document
print(doc, target = file.path(crowns_dir,"TABLE_S4_table_best_models_RATE_CHANGE_using_MEDIAN_VALUES.docx"))


## FIG 8:PLOT RATES OF CHANGE BASED ON FITTED MODELS OVER MEDIAN VALUES

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

best_models_median2 <- read.table(file.path(crowns_dir, "RATES_CHANGES_BEST_R2_MEDIAN_VALUES_all_sites.txt"), sep = "\t", header = TRUE)

#Define the desired order of responses
responses <-c("F1_cbh", "F1_depth", "F1_LAD","F1_Hdepth","MXLAD_cbh", "BP_cbh", "BR_cbh",
              "MXLAD_depth", "BP_depth", "MXLAD_dist", "BP_dist",
              "MXLAD_Hdepth","BP_Hdepth", "MAX_height", "total_LAI", "understory_LAI")

# Reorder factor levels for the Variable column based on the desired order
best_models_median2$Variable <- factor(best_models_median2$Variable, levels = responses)

# Reorder the rows based on 'Variable'
best_models_median2 <- best_models_median2[order(best_models_median2$Variable), ]

########  PLOTTING  ###############3

# Initialize a list to hold all individual plots across all groups
all_plots <- list()

# Set the output directory for saving plots
output_dir <- "DATA_AVAILABILITY/2_THIN_LOWEST/RESULTS/MODELS"

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
                          lower_ci = ~ quantile(.x, 0.25, na.rm = TRUE),
                          upper_ci = ~ quantile(.x, 0.75, na.rm = TRUE)),
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
  rsq_values <- best_models_median2 %>%
    dplyr::filter(Variable == var) %>%
    dplyr::select(Variable, rate_of_change, R_squared) %>%
    distinct() %>%
    mutate(rate_of_change = as.numeric(rate_of_change),
           R_squared = as.numeric(R_squared))
  
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

print(faceted_plot)

# Save the Faceted Plot
ggsave(
  filename = file.path(output_dir, "FIG8_PLOT_RATES_CHANGE_BEST_FITTED_MODELS_MEDIAN_VARIABLES_p2575.png"),
  plot = faceted_plot,
  width = 12,
  height = 8,
  dpi = 300
)


