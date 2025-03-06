# caitlinch/treelikeness_metrics/code/func_plotting.R
# Caitlin Cherryh 2025

# This file contains functions for data formatting and plotting

#### Experiment 1 data processing and plotting ####
exp1.metric.statistics.wrapper <- function(metric_name, df) {
  ## Function to take a test statistic and return, for each set of parameters,
  # summary statistics (median, mean, range, upper and lower quantiles)

  # Trim df to only that metric
  metric_df <- df[, c("row_id",
                      "uid",
                      "num_taxa",
                      "num_trees",
                      "tree_depth",
                      metric_name)]
  # Identify unique combinations of input parameters
  params_df <- unique(df[, c("num_taxa", "num_trees", "tree_depth")])
  # Iterate through rows of the params_df to extract summary statistics
  summary_list <- lapply(
    1:nrow(params_df),
    exp1.extract.summary.statistics,
    params_df = params_df,
    metric_df = metric_df,
    metric_name = metric_name
  )
  summary_df <- as.data.frame(do.call(rbind, summary_list))
  return(summary_df)
}



exp1.extract.summary.statistics <- function(
    row_num,
    params_df,
    metric_df,
    metric_name)
{
  ## Extract summary statistics for a single set of parameters
  params_row <- params_df[row_num, ]
  trimmed_df <- metric_df[which(
    metric_df$num_taxa == params_row$num_taxa &
      metric_df$num_trees == params_row$num_trees &
      metric_df$tree_depth == params_row$tree_depth
  ), ]
  result_vec <- trimmed_df[, metric_name]
  if (length(which(is.na(result_vec))) == length(result_vec)) {
    # If all values are NA, cannot calculate summary statistics
    output_row <- c(
      params_row$num_taxa[1],
      params_row$num_trees[1],
      params_row$tree_depth[1],
      metric_name,
      NA,
      NA,
      NA,
      NA,
      NA,
      NA,
      length(result_vec),
      length(which(is.na(result_vec))),
      TRUE
    )
  } else {
    output_row <- c(
      params_row$num_taxa[1],
      params_row$num_trees[1],
      params_row$tree_depth[1],
      metric_name,
      mean(result_vec, na.rm = TRUE),
      median(result_vec, na.rm = TRUE),
      min(result_vec, na.rm = TRUE),
      max(result_vec, na.rm = TRUE),
      summary(result_vec)[["1st Qu."]],
      summary(result_vec)[["3rd Qu."]],
      length(result_vec),
      length(which(is.na(result_vec))),
      FALSE
    )
  }
  names(output_row) <- c(
    "num_taxa",
    "num_trees",
    "tree_depth",
    "metric_name",
    "mean",
    "median",
    "minimum",
    "maximum",
    "first_quantile",
    "third_quantile",
    "num_values",
    "num_NA_values",
    "all_values_NA"
  )
  return(output_row)
}



exp1.process.NetworkTreelikenessTest <- function(exp1_df) {
  ## Transform the Network Treelikeness Test results into more plottable format
  # Make a table of all possible parameter values for the network treelikeness test
  ntlt_df <- expand.grid(
    "num_taxa" = unique(exp1_df$num_taxa),
    "num_trees" = unique(exp1_df$num_trees),
    "tree_depth" = unique(exp1_df$tree_depth)
  )
  # Calculate proportion of treelike alignments for each set of parameter values
  prop_tl_results <- unlist(
    lapply(
      1:nrow(ntlt_df),
      reformat.network.treelikeness.test.results.exp1,
      params_df = ntlt_df,
      results_df = exp1_df
    )
  )
  # Add columns to match the exp1_long_df
  op_df <- data.frame(
    "num_taxa" = ntlt_df$num_taxa,
    "num_trees" = ntlt_df$num_trees,
    "tree_depth" = ntlt_df$tree_depth,
    "metric_name" = "NetworkTreelikenessTest",
    "mean" = prop_tl_results,
    "median" = prop_tl_results,
    "minimum" = NA,
    "maximum" = NA,
    "first_quantile" = NA,
    "third_quantile" = NA,
    "num_values" = unlist(lapply(1:nrow(ntlt_df), function(i) {
      length(
        which(
          exp1_df$num_taxa == ntlt_df[i, ]$num_taxa &
            exp1_df$num_trees == ntlt_df[i, ]$num_trees &
            exp1_df$tree_depth == ntlt_df[i, ]$tree_depth
        )
      )
    })),
    "num_NA_values" = unlist(lapply(1:nrow(ntlt_df), function(i) {
      length(which(is.na(exp1_df[which(
        exp1_df$num_taxa == ntlt_df[i, ]$num_taxa &
          exp1_df$num_trees == ntlt_df[i, ]$num_trees &
          exp1_df$tree_depth == ntlt_df[i, ]$tree_depth
      ), "NetworkTreelikenessTest"])))
    }))
  )
  op_df$all_values_NA <- op_df$num_values == op_df$num_NA_values
  return(op_df)
}



exp1.format.summary.dataframe<- function(df){
  ## Reformat summary data frame for pretty plots
  # Convert columns to numeric
  i <- c(1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12)
  df[, i] <- apply(df[, i], 2, function(x) as.numeric(as.character(x)))
  # Convert column to logical
  df$all_values_NA <- as.logical(df$all_values_NA)
  # Format Number of Taxa as a factor
  df$num_taxa <- factor(df$num_taxa,
                        levels = c(5, 10, 20, 50, 100),
                        labels = c(5, 10, 20, 50, 100),
                        ordered = T)
  # Format Tree Depth as a factor
  df$tree_depth <- factor(df$tree_depth,
                          levels = c(0.01, 0.10, 1.00),
                          labels = c(0.01, 0.10, 1.00),
                          ordered = T)
  # Create factor column of pretty print labels with the "metric_name" column
  df$var_label <-
    factor(
      df$metric_name,
      levels = c(
        "tree_proportion",
        "Cunningham_test",
        "mean_TIGER_value",
        "mean_delta_plot_value",
        "mean_Q_residual",
        "sCFL_mean",
        "LM_proportion_resolved_quartets",
        "NetworkTreelikenessTest"
      ),
      ordered = TRUE,
      labels = c(
        expression(atop("Tree", "proportion")),
        expression(atop("Cunningham", "metric")),
        expression(atop("Mean", "TIGER value")),
        expression(paste('Mean ', delta["q"])),
        expression(atop("Mean", "Q-Residual value")),
        expression(atop("Mean", "sCF value")),
        expression(atop("Proportion", "resolved quartets")),
        expression(atop("Proportion", "treelike alignments"))
      )
    )
  return(df)
}



exp1.format.dataframe<- function(df){
  ## Reformat data frame for pretty plots
  # Remove unneeded columns
  df <- df[, c(
    "num_taxa",
    "num_trees",
    "tree_depth",
    "tree_proportion",
    "Cunningham_test",
    "mean_delta_plot_value",
    "LM_proportion_resolved_quartets",
    "mean_Q_residual",
    "sCFL_mean",
    "mean_TIGER_value",
    "NetworkTreelikenessTest"
  )]
  # Melt df without Network Treelikeness Test
  nonbinary_metric_df <- df[, c(
    "num_taxa",
    "num_trees",
    "tree_depth",
    "tree_proportion",
    "Cunningham_test",
    "mean_delta_plot_value",
    "LM_proportion_resolved_quartets",
    "mean_Q_residual",
    "sCFL_mean",
    "mean_TIGER_value"
  )]
  long_df <- melt(nonbinary_metric_df,
                  id.vars = c("num_taxa", "num_trees", "tree_depth"))
  names(long_df) <- c("num_taxa", "num_trees", "tree_depth", "metric_name", "value")
  # Extract NTLT df
  ntlt_df <- exp1.process.NetworkTreelikenessTest(exp1_df)[, c(1:4,6)]
  names(ntlt_df) <- c("num_taxa", "num_trees", "tree_depth", "metric_name", "value")
  # Bind data frames
  output_df <- rbind(long_df, ntlt_df)
  # Format Number of Taxa as a factor
  output_df$num_taxa <- factor(output_df$num_taxa,
                               levels = c(5, 10, 20, 50, 100),
                               labels = c(5, 10, 20, 50, 100),
                               ordered = T)
  # Format Tree Depth as a factor
  output_df$tree_depth <- factor(output_df$tree_depth,
                                 levels = c(0.01, 0.10, 1.00),
                                 labels = c(0.01, 0.10, 1.00),
                                 ordered = T)
  # Create factor column of variables
  output_df$metric_name <-
    factor(
      output_df$metric_name,
      levels = c(
        "tree_proportion",
        "Cunningham_test",
        "mean_TIGER_value",
        "mean_delta_plot_value",
        "mean_Q_residual",
        "sCFL_mean",
        "LM_proportion_resolved_quartets",
        "NetworkTreelikenessTest"
      ),
      ordered = TRUE,
      labels = c(
        "tree_proportion",
        "Cunningham_test",
        "mean_TIGER_value",
        "mean_delta_plot_value",
        "mean_Q_residual",
        "sCFL_mean",
        "LM_proportion_resolved_quartets",
        "NetworkTreelikenessTest"
      )
    )
  # Create factor column of pretty print labels with the "metric_name" column
  output_df$var_label <-
    factor(
      output_df$metric_name,
      levels = c(
        "tree_proportion",
        "Cunningham_test",
        "mean_TIGER_value",
        "mean_delta_plot_value",
        "mean_Q_residual",
        "sCFL_mean",
        "LM_proportion_resolved_quartets",
        "NetworkTreelikenessTest"
      ),
      ordered = TRUE,
      labels = c(
        expression(atop("Tree", "proportion")),
        expression(atop("Cunningham", "metric")),
        expression(atop("Mean", "TIGER value")),
        expression(paste('Mean ', delta["q"])),
        expression(atop("Mean", "Q-Residual value")),
        expression(atop("Mean", "sCF value")),
        expression(atop("Proportion", "resolved quartets")),
        expression(atop("Proportion", "treelike alignments"))
      )
    )
  return(output_df)
}


#### Experiment 3 data processing ####
exp3.metric.statistics.wrapper <- function(metric_name, df) {
  ## Function to take a test statistic and return, for each set of parameters,
  # summary statistics (median, mean, range, upper and lower quantiles)

  # Trim df to only that metric
  metric_df <- df[, c(
    "row_id",
    "uid",
    "num_taxa",
    "recombination_value",
    "recombination_type",
    "speciation_rate",
    "num_trees",
    "tree_depth_subspersite",
    "tree_age",
    metric_name
  )]
  # Identify unique combinations of input parameters
  params_df <- unique(df[, c(
    "num_taxa",
    "recombination_value",
    "recombination_type",
    "speciation_rate",
    "num_trees",
    "tree_depth_subspersite",
    "tree_age"
  )])
  # Iterate through rows of the params_df to extract summary statistics
  summary_list <- lapply(
    1:nrow(params_df),
    exp3.extract.summary.statistics,
    params_df = params_df,
    metric_df = metric_df,
    metric_name = metric_name
  )
  summary_df <- as.data.frame(do.call(rbind, summary_list))
  return(summary_df)
}



exp3.extract.summary.statistics <- function(
    row_num,
    params_df,
    metric_df,
    metric_name)
{
  ## Extract summary statistics for a single set of parameters
  params_row <- params_df[row_num, ]
  trimmed_df <- metric_df[which(
    metric_df$num_taxa == params_row$num_taxa &
      metric_df$recombination_value == params_row$recombination_value &
      metric_df$recombination_type == params_row$recombination_type &
      metric_df$speciation_rate == params_row$speciation_rate &
      metric_df$num_trees == params_row$num_trees &
      metric_df$tree_depth_subspersite == params_row$tree_depth_subspersite &
      metric_df$tree_age == params_row$tree_age
  ), ]
  result_vec <- trimmed_df[, metric_name]
  if (length(which(is.na(result_vec))) == length(result_vec)) {
    # If all values are NA, cannot calculate summary statistics
    output_row <- as.character(c(
      params_row,
      metric_name,
      NA,
      NA,
      NA,
      NA,
      NA,
      NA,
      length(result_vec),
      length(which(is.na(result_vec))),
      TRUE
    ))
  } else {
    output_row <- as.character(c(
      params_row,
      metric_name,
      mean(result_vec, na.rm = TRUE),
      median(result_vec, na.rm = TRUE),
      min(result_vec, na.rm = TRUE),
      max(result_vec, na.rm = TRUE),
      summary(result_vec)[["1st Qu."]],
      summary(result_vec)[["3rd Qu."]],
      length(result_vec),
      length(which(is.na(result_vec))),
      FALSE
    ))
  }
  names(output_row) <- c(
    names(params_row),
    "metric_name",
    "mean",
    "median",
    "minimum",
    "maximum",
    "first_quantile",
    "third_quantile",
    "num_values",
    "num_NA_values",
    "all_values_NA"
  )
  return(output_row)
}



exp3.process.NetworkTreelikenessTest <- function(exp3_df) {
  ## Transform the Network Treelikeness Test results into more plottable format
  # Make a table of all possible parameter values for the network treelikeness test
  ntlt_params <- unique(exp3_df[, c(
    "num_taxa",
    "tree_age",
    "speciation_rate",
    "tree_depth_subspersite",
    "recombination_value",
    "recombination_type"
  )])
  # Calculate proportion of treelike alignments for each set of parameter values
  prop_tl_results <- unlist(
    lapply(
      1:nrow(ntlt_params),
      reformat.network.treelikeness.test.results.exp3,
      params_df = ntlt_params,
      results_df = exp3_df
    )
  )
  # Format output
  op_df <- data.frame(
    "num_taxa" = ntlt_params$num_taxa,
    "recombination_value" = ntlt_params$recombination_value,
    "recombination_type" = ntlt_params$recombination_type,
    "speciation_rate" = ntlt_params$speciation_rate,
    "num_trees" = 200,
    "tree_depth_subspersite" = ntlt_params$tree_depth_subspersite,
    "tree_age" = ntlt_params$tree_age,
    "metric_name" = "NetworkTreelikenessTest",
    "mean" = prop_tl_results,
    "median" = prop_tl_results,
    "minimum" = NA,
    "maximum" = NA,
    "first_quantile" = NA,
    "third_quantile" = NA
  )
  # Check number of values per
  op_df$num_values = unlist(lapply(1:nrow(ntlt_params), function(i) {
    length(
      which(
        exp3_df$num_taxa == ntlt_params[i, ]$num_taxa &
          exp3_df$tree_age == ntlt_params[i, ]$tree_age &
          exp3_df$speciation_rate == ntlt_params[i, ]$speciation_rate &
          exp3_df$tree_depth_subspersite == ntlt_params[i, ]$tree_depth_subspersite &
          exp3_df$recombination_value == ntlt_params[i, ]$recombination_value &
          exp3_df$recombination_type == ntlt_params[i, ]$recombination_type
      )
    )
  }))
  op_df$num_NA_values = unlist(lapply(1:nrow(ntlt_params), function(i) {
    length(which(is.na(exp3_df[which(
      exp3_df$num_taxa == ntlt_params[i, ]$num_taxa &
        exp3_df$tree_age == ntlt_params[i, ]$tree_age &
        exp3_df$speciation_rate == ntlt_params[i, ]$speciation_rate &
        exp3_df$tree_depth_subspersite == ntlt_params[i, ]$tree_depth_subspersite &
        exp3_df$recombination_value == ntlt_params[i, ]$recombination_value &
        exp3_df$recombination_type == ntlt_params[i, ]$recombination_type
    ), "NetworkTreelikenessTest"])))
  }))
  op_df$all_values_NA <- op_df$num_values == op_df$num_NA_values
  return(op_df)
}


exp3.format.summary.dataframe<- function(df){
  ## Reformat summary data frame for pretty plots
  # Convert columns to numeric
  i <- c(1, 2, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16)
  df[, i] <- apply(df[, i], 2, function(x) as.numeric(as.character(x)))
  # Convert column to logical
  df$all_values_NA <- as.logical(df$all_values_NA)
  # Format Number of Taxa as a factor
  df$num_taxa <- factor(df$num_taxa,
                        levels = c(5, 10, 20, 50, 100),
                        labels = c(5, 10, 20, 50, 100),
                        ordered = T)
  # Format Tree Depth and Tree Age as a factor
  df$tree_depth_subspersite <- factor(df$tree_depth_subspersite,
                                      levels = c(0.005, 0.050, 0.500),
                                      labels = c(0.005, 0.050, 0.500),
                                      ordered = T)
  df$tree_age <- factor(df$tree_age,
                        levels = c(5, 50, 500),
                        labels = c(5, 50, 500),
                        ordered = T)
  # Create factor column of pretty print labels with the "metric_name" column
  df$var_label <-
    factor(
      df$metric_name,
      levels = c(
        "tree_proportion",
        "Cunningham_test",
        "mean_TIGER_value",
        "mean_delta_plot_value",
        "mean_Q_residual",
        "sCF_mean",
        "LM_proportion_resolved_quartets",
        "NetworkTreelikenessTest"
      ),
      ordered = TRUE,
      labels = c(
        expression(atop("Tree", "proportion")),
        expression(atop("Cunningham", "metric")),
        expression(atop("Mean", "TIGER value")),
        expression(paste('Mean ', delta["q"])),
        expression(atop("Mean", "Q-Residual value")),
        expression(atop("Mean", "sCF value")),
        expression(atop("Proportion", "resolved quartets")),
        expression(atop("Proportion", "treelike alignments"))
      )
    )
  return(df)
}



exp3.format.dataframe<- function(df){
  ## Reformat data frame for pretty plots
  # Remove unneeded columns
  df <- df[, c(
    "num_taxa",
    "recombination_value",
    "recombination_type",
    "speciation_rate",
    "num_trees",
    "tree_depth_subspersite",
    "tree_age",
    "tree_proportion",
    "Cunningham_test",
    "mean_delta_plot_value",
    "LM_proportion_resolved_quartets",
    "mean_Q_residual",
    "sCF_mean",
    "mean_TIGER_value",
    "NetworkTreelikenessTest"
  )]
  # Melt df without Network Treelikeness Test
  nonbinary_metric_df <- df[, c(
    "num_taxa",
    "recombination_value",
    "recombination_type",
    "speciation_rate",
    "num_trees",
    "tree_depth_subspersite",
    "tree_age",
    "tree_proportion",
    "Cunningham_test",
    "mean_delta_plot_value",
    "LM_proportion_resolved_quartets",
    "mean_Q_residual",
    "sCF_mean",
    "mean_TIGER_value"
  )]
  long_df <- melt(
    nonbinary_metric_df,
    id.vars = c(
      "num_taxa",
      "recombination_value",
      "recombination_type",
      "speciation_rate",
      "num_trees",
      "tree_depth_subspersite",
      "tree_age"
    )
  )
  names(long_df) <- c(
    "num_taxa",
    "recombination_value",
    "recombination_type",
    "speciation_rate",
    "num_trees",
    "tree_depth_subspersite",
    "tree_age",
    "metric_name",
    "value"
  )
  # Extract NTLT df
  ntlt_df <- exp3.process.NetworkTreelikenessTest(exp3_df)[, c(1:9)]
  names(ntlt_df) <- c("num_taxa", "recombination_value", "recombination_type",
                      "speciation_rate", "num_trees", "tree_depth_subspersite",
                      "tree_age", "metric_name", "value")
  # Bind data frames
  output_df <- rbind(long_df, ntlt_df)
  # Format Number of Taxa as a factor
  output_df$num_taxa <- factor(output_df$num_taxa,
                        levels = c(5, 10, 20, 50, 100),
                        labels = c(5, 10, 20, 50, 100),
                        ordered = T)
  # Format Tree Depth and Tree Age as a factor
  output_df$tree_depth_subspersite <- factor(output_df$tree_depth_subspersite,
                                      levels = c(0.005, 0.050, 0.500),
                                      labels = c(0.005, 0.050, 0.500),
                                      ordered = T)
  output_df$tree_age <- factor(output_df$tree_age,
                        levels = c(5, 50, 500),
                        labels = c(5, 50, 500),
                        ordered = T)
  # Create factor column of pretty print labels with the "metric_name" column
  output_df$var_label <-
    factor(
      output_df$metric_name,
      levels = c(
        "tree_proportion",
        "Cunningham_test",
        "mean_TIGER_value",
        "mean_delta_plot_value",
        "mean_Q_residual",
        "sCF_mean",
        "LM_proportion_resolved_quartets",
        "NetworkTreelikenessTest"
      ),
      ordered = TRUE,
      labels = c(
        expression(atop("Tree", "proportion")),
        expression(atop("Cunningham", "metric")),
        expression(atop("Mean", "TIGER value")),
        expression(paste('Mean ', delta["q"])),
        expression(atop("Mean", "Q-Residual value")),
        expression(atop("Mean", "sCF value")),
        expression(atop("Proportion", "resolved quartets")),
        expression(atop("Proportion", "treelike alignments"))
      )
    )
  return(output_df)
}

