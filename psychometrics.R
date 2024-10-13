# full code for scale development (needs checking)
# this shows a three factor solution to a procrastination scale,
# with item analysis code, exploratory factor analysis code, reliability measures, 
# validity measures (using concurrent, discriminatory validity), 
# and confirmatory factor analysis.



# packages ----------------------------------------------------------------

library(psych)
library(dplyr)
library(ggplot2)
library(reshape2)
library(lavaan)
library(semPlot)
library(tidyverse)
library(reshape2)
library(broom)

# Set options to avoid scientific notation
options(scipen = 999)



# Data preparation and transformation -------------------------------------


# List of items to reverse
reverse_items <- c(2, 6, 7, 9, 12, 13, 14, 15, 19, 22, 25, 29)

# Reverse the items, rename them as 'proc_X_reversed', and remove the original items
# the warning message will shown you which items have been deleted,
# these should be the original items before they have been reversed
df <- df_procrastination %>%
  # First, create the reversed items and rename them with "_reversed" suffix
  mutate_at(vars(paste0("proc_", reverse_items)), ~ 5 - .) %>%
  rename_at(vars(paste0("proc_", reverse_items)), ~ paste0("proc_", reverse_items, "_reversed")) %>%
  # Then, remove the original unreversed items
  select(-one_of(paste0("proc_", reverse_items))) %>%
  # Explicitly remove the column with the uppercase 'Proc_1_reversed'
  select(-Proc_1_reversed) #for some reason it was making a new variable with a capital

# View the updated dataframe
view(df)



# item analysis -----------------------------------------------------------


# Define a function for item analysis
item_analysis <- function(data) {
  # Get item statistics for all items
  item_stats <- describe(data)
  
  # Get the names of proc items
  proc_items <- names(data)[startsWith(names(data), "proc")]
  
  # Initialize item-rest correlations for proc items
  item_rest_correlations <- sapply(proc_items, function(item_name) {
    item_vector <- data[[item_name]]
    rest_vector <- rowMeans(data[, proc_items[proc_items != item_name]], na.rm = TRUE)
    cor(item_vector, rest_vector, use = "complete.obs")
  })
  
  # Calculate frequency of endorsement for all items
  endorsement_results <- sapply(names(data), function(item_name) {
    item_vector <- data[[item_name]]
    freq_table <- table(item_vector, useNA = "no") / length(item_vector) * 100  # Get percentages
    
    # Check for endorsement criterion: adjacent scale points should each have >= 20%
    endorsement_check <- FALSE
    if (length(freq_table) >= 5) {  # Ensure there are enough points for 5-point scale
      endorsement_check <- any(sapply(1:4, function(i) {
        return(freq_table[i] >= 20 && freq_table[i + 1] >= 20)
      }))
    }
    return(endorsement_check)
  })
  
  # Combine statistics for proc items
  item_analysis_results <- data.frame(
    Item = proc_items,
    Mean = item_stats$mean[proc_items],
    SD = item_stats$sd[proc_items],
    Min = item_stats$min[proc_items],
    Max = item_stats$max[proc_items],
    Item_Rest_Corr = item_rest_correlations,
    Endorsement_Meets_Criterion = endorsement_results[proc_items]  # Only include for proc items
  )
  
  # Combine statistics for all items
  all_items_results <- data.frame(
    Item = names(data),
    Mean = item_stats$mean,
    SD = item_stats$sd,
    Min = item_stats$min,
    Max = item_stats$max,
    Endorsement_Meets_Criterion = endorsement_results
  )
  
  return(list(
    Proc_Items_Statistics = item_analysis_results,
    All_Items_Statistics = all_items_results
  ))
}

# Run item analysis
item_analysis_results <- item_analysis(df)

# Print the results
print(item_analysis_results$Proc_Items_Statistics)
print(item_analysis_results$All_Items_Statistics)


# Remove specified items from the dataframe
df <- df %>% 
  select(-proc_3, -proc_4, -proc_8, -proc_11, -proc_20, -proc_12_reversed, 
         -proc_29_reversed, -proc_13_reversed)

# View the updated dataframe
view(df)

# EFA - Exploratory Factor Analysis ------------------------------------------------------
# Subset the dataframe to only include items that begin with "proc"
proc_items_df <- df %>% select(starts_with("proc"))
# View the selected proc items
view(proc_items_df)


# Check KMO
KMO(proc_items_df)

# Bartlett's Test of Sphericity
cortest.bartlett(proc_items_df)

#remove items based on KMO and Bartlett's Test:

proc_items_df <- proc_items_df %>% 
  select(-proc_19_reversed)

# Conduct EFA without specifying number of factors to explore eigenvalues first
efa_results <- fa(proc_items_df, nfactors = length(proc_items_df), rotate = "none", fm = "ml") 

# Extract eigenvalues
eigenvalues <- efa_results$values

# Create a data frame for eigenvalues
eigenvalues_df <- data.frame(Factor = 1:length(eigenvalues), Eigenvalue = eigenvalues)

# Count factors with eigenvalues > 1
num_factors_greater_than_1 <- sum(eigenvalues > 1)
cat("Number of factors with eigenvalues > 1:", num_factors_greater_than_1, "\n")

# Show eigenvalues greater than 1
eigenvalues_greater_than_1 <- eigenvalues[eigenvalues > 1]
cat("Eigenvalues greater than 1:\n")
print(eigenvalues_greater_than_1)

# Plot eigenvalues with a red line at 1
ggplot(eigenvalues_df, aes(x = Factor, y = Eigenvalue)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed") +  # Add red dashed line at y = 1
  labs(title = "Scree Plot", x = "Factors", y = "Eigenvalues") +
  theme_minimal() +
  annotate("text", x = max(eigenvalues_df$Factor), y = 1.5, label = "Eigenvalue = 1", color = "red", hjust = 1)

# Now, perform EFA with appropriate number of factors based on eigenvalues > 1
efa_results_final <- fa(proc_items_df, nfactors = num_factors_greater_than_1, rotate = "varimax", fm = "pa") # principal axis

# Extract and show loadings only for "proc" items
loadings <- efa_results_final$loadings

# Filter and show only "proc" loadings
proc_loadings <- loadings[grepl("^proc", rownames(loadings)), , drop = FALSE]
print(proc_loadings)

# creating the factor solution --------------------------------------------

# create new subscales (factors) from items:
# some of these will have already been removed however


F1 <- df %>% 
  select(proc_1, proc_5, proc_18, proc_21)

F2 <- df %>% 
  select(proc_9_reversed, proc_10, proc_25_reversed)

F3 <- df %>% 
  select(proc_14_reversed, proc_15_reversed, proc_28)



#dont include the average score in the reliability analyses

# Compute Cronbach's alpha for the proc items included in the final model

psych::alpha(df[, c("proc_1", "proc_5", "proc_18", "proc_21",
                    "proc_9_reversed", "proc_10", "proc_25_reversed",
                    "proc_14_reversed", "proc_15_reversed", "proc_28")])


# Compute omega reliability for the same items
omega(df[, c("proc_1", "proc_5", "proc_18", "proc_21",
             "proc_9_reversed", "proc_10", "proc_25_reversed",
             "proc_14_reversed", "proc_15_reversed", "proc_28")])

# validity testing --------------------------------------------------------

#create an avg score for each scale, and the new scale
df <- df %>%
  mutate(proc = rowMeans(select(., 
                                "proc_1", 
                                "proc_5", 
                                "proc_18", 
                                "proc_21", 
                                "proc_9_reversed", 
                                "proc_10", 
                                "proc_25_reversed", 
                                "proc_14_reversed", 
                                "proc_15_reversed", 
                                "proc_28"), 
                         na.rm = TRUE)) %>%
  mutate(cproc = rowMeans(select(., starts_with("cproc_")), na.rm = TRUE)) %>%
  mutate(dep = rowMeans(select(., starts_with("dep_")), na.rm = TRUE)) %>%
  mutate(sd = rowMeans(select(., starts_with("sd_")), na.rm = TRUE)) %>%
  mutate(mind = rowMeans(select(., starts_with("mind_")), na.rm = TRUE))


view(df)

# Create a data frame with the new variables
new_vars <- df %>% select(proc, cproc, dep, sd, mind)  # Select only the new variables

# Inspect the structure and summary statistics of new_vars
print("Structure of new_vars:")
str(new_vars)
print("Summary statistics of new_vars:")
summary(new_vars)

# Check for NA values in the new variables
na_count <- colSums(is.na(new_vars))
print("Count of NA values in each variable:")
print(na_count)

# Optionally remove rows with NA values in the new variables
cleaned_new_vars <- new_vars %>% drop_na()

# Create a matrix to hold p-values
p_matrix <- matrix(NA, ncol = ncol(cleaned_new_vars), nrow = ncol(cleaned_new_vars))  # Initialize p_matrix

# Compute p-values for each pair of variables
for (i in 1:ncol(cleaned_new_vars)) {
  for (j in 1:ncol(cleaned_new_vars)) {
    if (i != j) {
      # Perform correlation test and check if the variables have enough non-NA observations
      valid_obs <- complete.cases(cleaned_new_vars[[i]], cleaned_new_vars[[j]])  # Identify valid observations
      if (sum(valid_obs) > 2) {  # Ensure enough observations
        test <- cor.test(cleaned_new_vars[[i]][valid_obs], cleaned_new_vars[[j]][valid_obs])  # Perform correlation test
        p_matrix[i, j] <- test$p.value  # Store p-value
        # Format the p-value to 3 decimal places without scientific notation
        formatted_p_value <- format(p_matrix[i, j], nsmall = 3, scientific = FALSE)
        cat("Correlation between", names(cleaned_new_vars)[i], "and", names(cleaned_new_vars)[j], "- p-value:", formatted_p_value, "\n")
      } else {
        p_matrix[i, j] <- NA  # Not enough valid observations for correlation test
        cat("Not enough valid observations for", names(cleaned_new_vars)[i], "and", names(cleaned_new_vars)[j], "\n")
      }
    } else {
      p_matrix[i, j] <- NA  # No p-value for self-correlation
    }
  }
}

# Create the correlation matrix
cor_matrix <- cor(cleaned_new_vars, use = "complete.obs")  # Compute the correlation matrix

# Print the correlation matrix
print("Correlation Matrix:")
print(cor_matrix)

# Print the p-value matrix with formatted p-values
formatted_p_matrix <- format(p_matrix, nsmall = 3, scientific = FALSE)
print("P-value Matrix:")
print(formatted_p_matrix)


# confirmatory factor analysis --------------------------------------------

psych::describe(df)


# Specify the CFA model for proc items
cfa_model <- '
  F1 =~ proc_1 + proc_5 + proc_18 + proc_21
  F2 =~ proc_9_reversed + proc_10 + proc_25_reversed
  F3 =~ proc_14_reversed + proc_15_reversed + proc_28
  proc_CFA =~ F1 + F2 + F3
'


# Fit the CFA model using the dataframe
cfa_results <- cfa(cfa_model, data = df, estimator = "MLR") #for non-normal and missing data 

# Summarize the CFA results
summary(cfa_results, fit.measures = TRUE, standardized = TRUE)

# Visualize the CFA model
semPaths(cfa_results, what = "std", edge.label.cex = 0.8, fade = FALSE,
         layout = "tree", sizeMan = 5, sizeLat = 8, nCharNodes = 0,
         edge.color = "black", label.color = "black", label.cex = 0.8)