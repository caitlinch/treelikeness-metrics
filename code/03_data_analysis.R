# caitlinch/treelikeness-metrics/code/03_data_analysis.R
# Caitlin Cherryh 2025

# This program takes results from applying various treelikeness tests and
# performs data analysis/plotting

# Use the treelikeness-metrics RProject for this script



#### 1. Set parameters ####
# data_directory <- Directory containing result csv files
# plot_directory <- Directory for output of data analysis/plots
data_directory <- "output/"
plot_directory <- "plot/"



#### 2. Prepare analyses ####
# Source function files
source("code/func_data_analysis.R")

# Open packages
library(scales)
library(reshape2)
library(ggplot2)
library(patchwork)

# List all data files
data_files <- paste0(data_directory, list.files(data_directory))



#### 3. Prepare data from Experiment 1 for plotting ####
# Open data file from Experiment 1 as a dataframe
exp1_data_file <- grep(
  "exp1",
  grep("treelikeness_metrics_collated_results", data_files, value = TRUE),
  value = TRUE
)
exp1_df <- read.csv(file = exp1_data_file, stringsAsFactors = FALSE)

# Convert sCFL values to decimal from percentage
exp1_df$sCFL_mean <- exp1_df$sCFL_mean / 100

# Convert mean tiger value to numeric
# Note: 4% (150/3750) of simulations for experiment 1 do not have a TIGER value (TIGER failed to run in these cases)
#       Therefore converting to numeric will coerce these values (i.e. mean_TIGER_value = "no_TIGER_run") to NA
exp1_df$mean_TIGER_value <- as.numeric(exp1_df$mean_TIGER_value)

# Remove columns you don't want for plotting
nonbinary_metrics <- c(
  "tree_proportion",
  "Cunningham_test",
  "mean_delta_plot_value",
  "LM_proportion_resolved_quartets",
  "mean_Q_residual",
  "sCFL_mean",
  "mean_TIGER_value"
)
exp1_df <- exp1_df[, c(
  "row_id",
  "uid",
  "num_taxa",
  "num_trees",
  "tree_depth",
  nonbinary_metrics,
  "NetworkTreelikenessTest"
)]

# Create data frame of summary statistics for the non-binary metrics
exp1_summary_df <- as.data.frame(do.call(
  rbind,
  lapply(nonbinary_metrics, exp1.metric.statistics.wrapper, df = exp1_df)
))

# Reformat and bind the Network Treelikeness Test results
exp1_summary_df <- rbind(exp1_summary_df,
                         exp1.process.NetworkTreelikenessTest(exp1_df))

# Format the experiment 1 results for pretty plots
exp1_df <- exp1.format.dataframe(exp1_df)
exp1_summary_df <- exp1.format.summary.dataframe(exp1_summary_df)

# Replace any Cunningham metric values <0 with 0
# (for plotting purposes - to manage out of bounds values with free y axes)
exp1_summary_squish_df <- exp1_summary_df
exp1_summary_squish_df[which(exp1_summary_squish_df$minimum < 0), "minimum"] <- 0
exp1_summary_squish_df[which(exp1_summary_squish_df$maximum > 1), "maximum"] <- 1



#### 4. Plot data from Experiment 1 ####
# Set log10 minor breaks for x axis
x_axis_minor_breaks <-  unique(c(
  seq(1, 10, 1),
  seq(10, 100, 10),
  seq(100, 1000, 100),
  seq(1000, 10000, 1000)
))

## Plot 1:
# Lines showing median values for each test statistic as the number of trees
# increases, faceted by tree depth
#   - Log10 x-axis
#   - Free y-axis for each test statistic
#   - Use "squish" version of df to avoid the 2 minimum Cunningham values < 0
exp1_plot1 <-
  ggplot(exp1_summary_squish_df) +
  geom_ribbon(aes(
    x = num_trees,
    ymin = minimum,
    ymax = maximum,
    fill = as.factor(num_taxa)
  ),
  alpha = 0.2) +
  geom_line(aes(
    x = num_trees,
    y = median,
    color = as.factor(num_taxa)
  )) +
  facet_grid(var_label ~ tree_depth, scales = "free", labeller = label_parsed) +
  scale_x_log10(minor_breaks = x_axis_minor_breaks) +
  scale_y_continuous(name = "Median test statistic value", oob = scales::squish) +
  scale_color_viridis_d(direction = -1) +
  scale_fill_viridis_d(direction = -1) +
  guides(
    color = guide_legend(title = "Number of\ntaxa"),
    fill = guide_legend(title = "Number of\ntaxa")
  ) +
  labs(
    title = "Random Tree Simulations",
    subtitle = "Tree depth (substitutions per site)",
    x = expression("Number of trees (" * log[10] * " scale)")
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(
      size = 18,
      margin = margin(t = 15, r = 0, b = 0, l = 0)
    ),
    axis.text.x = element_text(
      size = 14,
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    axis.title.y = element_text(
      size = 18,
      margin = margin(t = 0, r = 15, b = 0, l = 0)
    ),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 10),
    plot.title = element_text(
      hjust = 0.5,
      size = 25,
      margin = margin(t = 0, r = 0, b = 15, l = 0)
    ),
    plot.subtitle = element_text(hjust = 0.5, size = 18)
  )
ggsave(exp1_plot1,
       filename = paste0(plot_directory, "exp1_plot1_median_ribbon.pdf"),
       width = 10,
       height = 13,
       units = "in"
)
# ggsave(exp1_plot1,
#        filename = paste0(plot_directory, "exp1_plot1_median_ribbon.png"),
#        width = 10,
#        height = 13,
#        units = "in"
# )


## Plot 2:
# Lines showing median values for each test statistic as the number of trees
# increases, faceted by tree depth
#   - Log10 x-axis
#   - Free y-axis for each test statistic and each tree depth
#   - Use "squish" version of df to avoid the 2 minimum Cunningham values < 0
exp1_plot2_panel1 <-
  ggplot(exp1_summary_squish_df[which(exp1_summary_squish_df$tree_depth == 0.01),]) +
  geom_ribbon(aes(
    x = num_trees,
    ymin = minimum,
    ymax = maximum,
    fill = as.factor(num_taxa)
  ),
  alpha = 0.2) +
  geom_line(aes(
    x = num_trees,
    y = median,
    color = as.factor(num_taxa)
  )) +
  facet_grid(var_label ~ tree_depth, scales = "free", labeller = label_parsed) +
  scale_x_log10(minor_breaks = x_axis_minor_breaks) +
  scale_y_continuous(name = "Median test statistic value", oob = scales::squish) +
  scale_color_viridis_d(direction = -1) +
  scale_fill_viridis_d(direction = -1) +
  guides(
    color = guide_legend(title = "Number of\ntaxa"),
    fill = guide_legend(title = "Number of\ntaxa")
  ) +
  labs(
    title = "Random Tree Simulations",
    subtitle = "Tree depth (substitutions per site)",
    x = expression("Number of trees (" * log[10] * " scale)")
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      size = 14,
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    axis.title.y = element_text(
      size = 18,
      margin = margin(t = 0, r = 15, b = 0, l = 0)
    ),
    axis.text.y = element_text(size = 14),
    legend.position = "none",
    strip.text.x = element_text(size = 11),
    strip.text.y = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

exp1_plot2_panel2 <-
  ggplot(exp1_summary_squish_df[which(exp1_summary_squish_df$tree_depth == 0.1),]) +
  geom_ribbon(aes(
    x = num_trees,
    ymin = minimum,
    ymax = maximum,
    fill = as.factor(num_taxa)
  ),
  alpha = 0.2) +
  geom_line(aes(
    x = num_trees,
    y = median,
    color = as.factor(num_taxa)
  )) +
  facet_grid(var_label ~ tree_depth, scales = "free", labeller = label_parsed) +
  scale_x_log10(minor_breaks = x_axis_minor_breaks) +
  scale_y_continuous(name = "Median test statistic value", oob = scales::squish) +
  scale_color_viridis_d(direction = -1) +
  scale_fill_viridis_d(direction = -1) +
  guides(
    color = guide_legend(title = "Number of\ntaxa"),
    fill = guide_legend(title = "Number of\ntaxa")
  ) +
  labs(
    title = "Random Tree Simulations",
    subtitle = "Tree depth (substitutions per site)",
    x = expression("Number of trees (" * log[10] * " scale)")
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(
      size = 18,
      margin = margin(t = 15, r = 0, b = 0, l = 0)
    ),
    axis.text.x = element_text(
      size = 14,
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 14),
    legend.position = "none",
    strip.text.x = element_text(size = 11),
    strip.text.y = element_blank(),
    plot.title = element_text(
      hjust = 0.5,
      size = 25,
      margin = margin(t = 0, r = 0, b = 15, l = 0)
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 18,
      margin = margin(t = 0, r = 0, b = 15, l = 0)
    )
  )

exp1_plot2_panel3 <-
  ggplot(exp1_summary_squish_df[which(exp1_summary_squish_df$tree_depth == 1),]) +
  geom_ribbon(aes(
    x = num_trees,
    ymin = minimum,
    ymax = maximum,
    fill = as.factor(num_taxa)
  ),
  alpha = 0.2) +
  geom_line(aes(
    x = num_trees,
    y = median,
    color = as.factor(num_taxa)
  )) +
  facet_grid(var_label ~ tree_depth, scales = "free", labeller = label_parsed) +
  scale_x_log10(minor_breaks = x_axis_minor_breaks) +
  scale_y_continuous(name = "Median test statistic value", oob = scales::squish) +
  scale_color_viridis_d(direction = -1) +
  scale_fill_viridis_d(direction = -1) +
  guides(
    color = guide_legend(title = "Number of\ntaxa"),
    fill = guide_legend(title = "Number of\ntaxa")
  ) +
  labs(
    title = "Random Tree Simulations",
    subtitle = "Tree depth (substitutions per site)",
    x = expression("Number of trees (" * log[10] * " scale)")
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      size = 14,
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    strip.text.x = element_text(size = 11),
    strip.text.y = element_text(size = 11),
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

exp1_plot2 <-
  exp1_plot2_panel1 + exp1_plot2_panel2 + exp1_plot2_panel3
ggsave(exp1_plot2,
       filename = paste0(plot_directory, "exp1_plot2_median_ribbon_freeY.pdf"),
       width = 12,
       height = 14,
       units = "in"
)
# ggsave(exp1_plot2,
#        filename = paste0(plot_directory, "exp1_plot2_median_ribbon_freeY.png"),
#        width = 12,
#        height = 14,
#        units = "in"
# )

## Plot 3:
# Dot plot of raw data
#   - Log10 x-axis
#   - Free y-axis for each column
# Note: 150 NAs from TIGER, 2 values < 0 from Cunningham metric
exp1_plot3 <-
  ggplot(exp1_df,
         aes(x = num_trees, y = value, color = as.factor(num_taxa)) ) +
  geom_point(alpha = 0.4) +
  facet_grid(var_label ~ tree_depth, scales = "free", labeller = label_parsed) +
  scale_x_log10(minor_breaks = x_axis_minor_breaks) +
  scale_y_continuous(
    name = "Test statistic value",
    limits = c(0,1.10),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = c(0, 0.25, 0.5, 0.75, 1),
    oob = scales::rescale_none) +
  scale_color_viridis_d(direction = -1) +
  guides(color = guide_legend(title = "Number of\ntaxa")) +
  labs(title = "Tree depth\n(substitutions per site)",
       x = expression("Number of trees (" * log[10] * " scale)")) +
  theme_bw() +
  theme(
    axis.title.x = element_text(
      size = 18,
      margin = margin(
        t = 15,
        r = 0,
        b = 0,
        l = 0
      )
    ),
    axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
    axis.title.y = element_text(
      size = 18,
      margin = margin(
        t = 0,
        r = 15,
        b = 0,
        l = 0
      )
    ),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 18)
  )
ggsave(exp1_plot3,
       filename = paste0(plot_directory, "exp1_plot3_points.pdf"),
       width = 10,
       height = 13,
       units = "in"
)
# ggsave(exp1_plot3,
#        filename = paste0(plot_directory, "exp1_plot3_points.png"),
#        width = 10,
#        height = 13,
#        units = "in"
# )


#### 5. Prepare data from Experiment 3 for plotting ####
# Open data file from Experiment 1 as a dataframe
exp3_data_file <- grep("00_",
                       grep(
                         "exp3",
                         grep("treelikeness_metrics_collated_results", data_files, value = TRUE),
                         value = TRUE
                       ),
                       value = TRUE,
                       invert = TRUE)
exp3_df <- read.csv(file = exp3_data_file, stringsAsFactors = FALSE)

# Convert sCFL values to decimal from percentage
exp1_df$sCFL_mean <- exp1_df$sCFL_mean / 100

# Convert mean tiger value to numeric
# Note: 4% (150/3750) of simulations for experiment 1 do not have a TIGER value (TIGER failed to run in these cases)
#       Therefore converting to numeric will coerce these values (i.e. mean_TIGER_value = "no_TIGER_run") to NA
exp1_df$mean_TIGER_value <- as.numeric(exp1_df$mean_TIGER_value)

# Remove columns you don't want for plotting
nonbinary_metrics <- c(
  "tree_proportion",
  "Cunningham_test",
  "mean_delta_plot_value",
  "LM_proportion_resolved_quartets",
  "mean_Q_residual",
  "sCFL_mean",
  "mean_TIGER_value"
)
exp1_df <- exp1_df[, c(
  "row_id",
  "uid",
  "num_taxa",
  "num_trees",
  "tree_depth",
  nonbinary_metrics,
  "NetworkTreelikenessTest"
)]

# Create data frame of summary statistics for the non-binary metrics
exp1_summary_df <- as.data.frame(do.call(
  rbind,
  lapply(nonbinary_metrics, exp1.metric.statistics.wrapper, df = exp1_df)
))

# Reformat and bind the Network Treelikeness Test results
exp1_summary_df <- rbind(exp1_summary_df,
                         exp1.process.NetworkTreelikenessTest(exp1_df))

# Format the experiment 1 results for pretty plots
exp1_df <- exp1.format.dataframe(exp1_df)
exp1_summary_df <- exp1.format.summary.dataframe(exp1_summary_df)

# Replace any Cunningham metric values <0 with 0
# (for plotting purposes - to manage out of bounds values with free y axes)
exp1_summary_squish_df <- exp1_summary_df
exp1_summary_squish_df[which(exp1_summary_squish_df$minimum < 0), "minimum"] <- 0
exp1_summary_squish_df[which(exp1_summary_squish_df$maximum > 1), "maximum"] <- 1



if (plot_exp3 == TRUE){
  ## Open data file from Experiment 3 as a dataframe
  exp3_data_file <- grep("00_", grep("exp3", grep("treelikeness_metrics_collated_results", data_files, value = TRUE), value = TRUE), value = TRUE, invert = TRUE)
  exp3_df <- read.csv(file = exp3_data_file, stringsAsFactors = FALSE)

  ## Manage test statistic columns
  # Convert sCFL values to decimal from percentage
  exp3_df$sCFL_mean <- exp3_df$sCF_mean / 100
  exp3_df$sCFL_median <- exp3_df$sCF_median / 100

  ## Melt exp3_wide_df for better plotting
  id_var_columns <- c("row_id", "uid", "num_trees", "num_taxa", "tree_age", "recombination_value", "recombination_type", "speciation_rate", "tree_depth_subspersite")
  exp3_long_df <- melt(exp3_df, id.vars = id_var_columns, measure.vars = c("LM_proportion_resolved_quartets", "sCFL_mean", "mean_delta_plot_value", "mean_Q_residual",
                                                                           "mean_TIGER_value", "Cunningham_test", "tree_proportion") )

  # Transform the Network Treelikeness Test results into more plottable format
  # Make a table of all possible parameter values for the network treelikeness test
  ntlt_params <- expand.grid("num_taxa" = sort(unique(exp3_df$num_taxa)), "tree_age" = sort(unique(exp3_df$tree_age)),
                             "speciation_rate" = sort(unique(exp3_df$speciation_rate)), "tree_depth_subspersite" = sort(unique(exp3_df$tree_depth_subspersite)),
                             "recombination_value" = sort(unique(exp3_df$recombination_value)), "recombination_type" = unique(exp3_df$recombination_type))
  # Calculate proportion of treelike alignments for each set of parameter values
  prop_tl_results <- unlist(lapply(1:nrow(ntlt_params), reformat.network.treelikeness.test.results.exp3, params_df = ntlt_params, results_df = exp3_df))
  # Add columns to match the exp3_long_df
  ntlt_params$row_id <- rep(NA, length(prop_tl_results))
  ntlt_params$uid <- rep(NA, length(prop_tl_results))
  ntlt_params$num_trees <- 1000
  ntlt_params$value <- prop_tl_results
  ntlt_params$variable <- "NetworkTreelikenessTest"
  # Restructure the dataframe to match the exp3_long_df
  ntlt_params <- ntlt_params[,c(names(exp3_long_df)),]
  # Remove rows with 5 taxa and Ancient event (not conducted)
  remove_rows <- which(ntlt_params$num_taxa == 5 & ntlt_params$recombination_type == "Ancient")
  keep_rows <- setdiff(1:nrow(ntlt_params), remove_rows)
  ntlt_params <- ntlt_params[keep_rows, ]
  # Bind to the exp3_long_df
  exp3_long_df <- rbind(exp3_long_df, ntlt_params)

  # Add fancy labels for facets (based on which test statistics were selected)
  exp3_long_df$var_label <- factor(exp3_long_df$variable,
                                   levels = c("tree_proportion", "Cunningham_test", "mean_delta_plot_value",
                                              "LM_proportion_resolved_quartets","NetworkTreelikenessTest",
                                              "mean_Q_residual", "sCFL_mean", "mean_TIGER_value"),
                                   ordered = TRUE,
                                   labels = c(expression(atop("Tree","proportion")), expression(atop("Cunningham","metric")),
                                              expression(paste('Mean ', delta["q"])), expression(atop("Proportion","resolved quartets")),
                                              expression(atop("Proportion","treelike alignments")), expression(atop("Mean", "Q-Residual value")),
                                              expression(atop("Mean", "sCF value")), expression(atop("Mean","TIGER value"))) )
} # end if (plot_exp3 == TRUE)



#### 6. Plot data from Experiment 3 ####
if (plot_exp3 == TRUE){
  # Set color palette for plots 1 and 2
  viridis_picks = scales::viridis_pal(direction = -1)(5)

  speciation_rates <- unique(exp3_long_df$speciation_rate)
  for (s in speciation_rates){
    print(s)
    # Remove speciation rate of 1.0
    s_df <- exp3_long_df[exp3_long_df$speciation_rate == s, ]

    ## Plot 1: Ancient events. Smooth lines showing average values for each test statistic as the amount of recombination increases, faceted by tree depth ##
    # Set dataset for plot
    plot_df <- s_df
    plot_df <- plot_df[plot_df$recombination_type == "Ancient", ]
    # Set plot title
    plot_title <- paste0("Ancient introgression; speciation rate = ", s)
    # Construct plot with fixed y axis
    p <- ggplot(plot_df, aes(x = recombination_value, y = value, color = as.factor(num_taxa))) +
      geom_smooth(method = "loess", alpha = 0.2, linewidth = 0, span = 0.75,
                  aes(ymin = ifelse(after_stat(ymin) < 0, 0, after_stat(ymin)), ymax = ifelse(after_stat(ymax) > 1, 1, after_stat(ymax)))) +
      stat_smooth(method = "loess", geom = "line", linewidth = 1.1, alpha = 0.7, span = 0.75) +
      facet_grid(var_label~tree_age, scales = "fixed", labeller = label_parsed) +
      scale_x_continuous(name = "Proportion of recombinant DNA", breaks = seq(0,0.5, 0.1), labels = seq(0,0.5, 0.1), minor_breaks = seq(0,0.5, 0.05)) +
      scale_y_continuous(name = "Test statistic value", limits = c(0,1.10), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1), oob=scales::rescale_none) +
      scale_color_manual(values = viridis_picks[2:5]) +
      guides(color = guide_legend(title = "Number of\ntaxa")) +
      labs(title = plot_title,
           subtitle = "Tree age (Million years)") +
      theme_bw() +
      theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
            axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 14),
            legend.title = element_text(size = 18, hjust = 0.5), legend.text = element_text(size = 16),
            strip.text = element_text(size = 11),
            plot.title = element_text(hjust = 0.5, size = 25, margin = margin(t = 0, r = 0, b = 15, l = 0)),
            plot.subtitle = element_text(size = 18, hjust = 0.5))
    # Save plot
    plot_prefix <- "mainFig_exp3_plot1_tree_depth_Ancient_s"
    ggsave(p, filename = paste0(plot_directory, "/pdf_plots/", plot_prefix, s, ".pdf"), width = 10, height = 12.5, units = "in")
    ggsave(p, filename = paste0(plot_directory, "/png_plots/", plot_prefix, s, ".png"), width = 10, height = 12.5, units = "in")

    # Construct plot with points only
    p <- ggplot(plot_df, aes(x = recombination_value, y = value, color = as.factor(num_taxa))) +
      geom_point(alpha = 0.4) +
      facet_grid(var_label~tree_age, scales = "fixed", labeller = label_parsed) +
      scale_x_continuous(name = "Proportion of recombinant DNA", breaks = seq(0,0.5, 0.1), labels = seq(0,0.5, 0.1), minor_breaks = seq(0,0.5, 0.05)) +
      scale_y_continuous(name = "Test statistic value", limits = c(0,1.10), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1), oob=scales::rescale_none) +
      scale_color_manual(values = viridis_picks[2:5]) +
      guides(color = guide_legend(title = "Number of\ntaxa")) +
      labs(title = "Tree age (Million years)") +
      theme_bw() +
      theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
            axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 14),
            legend.title = element_text(size = 18, hjust = 0.5), legend.text = element_text(size = 16),
            strip.text = element_text(size = 11),
            plot.title = element_text(size = 18, hjust = 0.5))
    # Save plot
    plot_prefix <- "exp3_plot1_tree_depth_Ancient_points_s"
    ggsave(p, filename = paste0(plot_directory, "/pdf_plots/", plot_prefix, s, ".pdf"), width = 10, height = 12.5, units = "in")
    ggsave(p, filename = paste0(plot_directory, "/png_plots/", plot_prefix, s, ".png"), width = 10, height = 12.5, units = "in")

    # Construct plot with free y axis
    p <- ggplot(plot_df, aes(x = recombination_value, y = value, color = as.factor(num_taxa))) +
      geom_smooth(method = "loess", alpha = 0.3, linewidth = 0, span = 0.75,
                  aes(ymin = ifelse(after_stat(ymin) < 0, 0, after_stat(ymin)), ymax = ifelse(after_stat(ymax) > 1, 1, after_stat(ymax)))) +
      stat_smooth(method = "loess", geom = "line", linewidth = 1.1, alpha = 1, span = 0.75) +
      facet_grid(var_label~tree_age, scales = "free_y", labeller = label_parsed) +
      scale_x_continuous(name = "Proportion of recombinant DNA", breaks = seq(0,0.5, 0.1), labels = seq(0,0.5, 0.1), minor_breaks = seq(0,0.5, 0.05)) +
      scale_y_continuous(name = "Test statistic value", oob=scales::rescale_none) +
      scale_color_manual(values = viridis_picks[2:5]) +
      guides(color = guide_legend(title = "Number of\ntaxa")) +
      labs(title = "Tree age (Million years)") +
      theme_bw() +
      theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
            axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 14),
            legend.title = element_text(size = 18, hjust = 0.5), legend.text = element_text(size = 16),
            strip.text = element_text(size = 11),
            plot.title = element_text(size = 18, hjust = 0.5))
    # Save plot
    plot_prefix <- "exp3_plot1_tree_depth_Ancient_freey_s"
    ggsave(p, filename = paste0(plot_directory, "/pdf_plots/", plot_prefix, s, ".pdf"), width = 10, height = 12.5, units = "in")
    ggsave(p, filename = paste0(plot_directory, "/png_plots/", plot_prefix, s, ".png"), width = 10, height = 12.5, units = "in")


    ## Plot 2: Recent events. Smooth lines showing average values for each test statistic as the amount of recombination increases, faceted by tree depth ##
    # Set dataset for plot
    plot_df <- s_df
    plot_df <- plot_df[plot_df$recombination_type == "Recent", ]
    # Set plot title
    plot_title <- paste0("Recent introgression; speciation rate = ", s)
    # Construct plot
    p <- ggplot(plot_df, aes(x = recombination_value, y = value, color = as.factor(num_taxa))) +
      geom_smooth(method = "loess", alpha = 0.2, linewidth = 0, span = 0.75,
                  aes(ymin = ifelse(after_stat(ymin) < 0, 0, after_stat(ymin)), ymax = ifelse(after_stat(ymax) > 1, 1, after_stat(ymax)))) +
      stat_smooth(method = "loess", geom = "line", linewidth = 1.1, alpha = 0.7, span = 0.75) +
      facet_grid(var_label~tree_age, scales = "fixed", labeller = label_parsed) +
      scale_x_continuous(name = "Proportion of recombinant DNA", breaks = seq(0,0.5, 0.1), labels = seq(0,0.5, 0.1), minor_breaks = seq(0,0.5, 0.05)) +
      scale_y_continuous(name = "Test statistic value", limits = c(0,1.10), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1), oob=scales::rescale_none) +
      scale_color_manual(values = viridis_picks[1:5]) +
      guides(color = guide_legend(title = "Number of\ntaxa")) +
      labs(title = plot_title,
           subtitle = "Tree age (Million years)") +
      theme_bw() +
      theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
            axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 14),
            legend.title = element_text(size = 18, hjust = 0.5), legend.text = element_text(size = 16),
            strip.text = element_text(size = 11),
            plot.title = element_text(hjust = 0.5, size = 25, margin = margin(t = 0, r = 0, b = 15, l = 0)),
            plot.subtitle = element_text(size = 18, hjust = 0.5))
    # Save plot
    plot_prefix <- "mainFig_exp3_plot2_tree_depth_Recent_s"
    ggsave(p, filename = paste0(plot_directory, "/pdf_plots/", plot_prefix, s, ".pdf"), width = 10, height = 12.5, units = "in")
    ggsave(p, filename = paste0(plot_directory, "/png_plots/", plot_prefix, s, ".png"), width = 10, height = 12.5, units = "in")

    # Construct plot with points only
    p <- ggplot(plot_df, aes(x = recombination_value, y = value, color = as.factor(num_taxa))) +
      geom_point(alpha = 0.4) +
      facet_grid(var_label~tree_age, scales = "fixed", labeller = label_parsed) +
      scale_x_continuous(name = "Proportion of recombinant DNA", breaks = seq(0,0.5, 0.1), labels = seq(0,0.5, 0.1), minor_breaks = seq(0,0.5, 0.05)) +
      scale_y_continuous(name = "Test statistic value", limits = c(0,1.10), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1), oob=scales::rescale_none) +
      scale_color_manual(values = viridis_picks[1:5]) +
      guides(color = guide_legend(title = "Number of\ntaxa")) +
      labs(title = "Tree age (Million years)") +
      theme_bw() +
      theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
            axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 14),
            legend.title = element_text(size = 18, hjust = 0.5), legend.text = element_text(size = 16),
            strip.text = element_text(size = 11),
            plot.title = element_text(size = 18, hjust = 0.5))
    # Save plot
    plot_prefix <- "exp3_plot2_tree_depth_Recent_points_s"
    ggsave(p, filename = paste0(plot_directory, "/pdf_plots/", plot_prefix, s, ".pdf"), width = 10, height = 12.5, units = "in")
    ggsave(p, filename = paste0(plot_directory, "/png_plots/", plot_prefix, s, ".png"), width = 10, height = 12.5, units = "in")

    # Construct plot with free y axis
    p <- ggplot(plot_df, aes(x = recombination_value, y = value, color = as.factor(num_taxa))) +
      geom_smooth(method = "loess", alpha = 0.3, linewidth = 0, span = 0.75,
                  aes(ymin = ifelse(after_stat(ymin) < 0, 0, after_stat(ymin)), ymax = ifelse(after_stat(ymax) > 1, 1, after_stat(ymax)))) +
      stat_smooth(method = "loess", geom = "line", linewidth = 1.1, alpha = 1, span = 0.75) +
      facet_grid(var_label~tree_age, scales = "free_y", labeller = label_parsed) +
      scale_x_continuous(name = "Proportion of recombinant DNA", breaks = seq(0,0.5, 0.1), labels = seq(0,0.5, 0.1), minor_breaks = seq(0,0.5, 0.05)) +
      scale_y_continuous(name = "Test statistic value", oob=scales::rescale_none) +
      scale_color_manual(values = viridis_picks[1:5]) +
      guides(color = guide_legend(title = "Number of\ntaxa")) +
      labs(title = "Tree age (Million years)") +
      theme_bw() +
      theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
            axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 14),
            legend.title = element_text(size = 18, hjust = 0.5), legend.text = element_text(size = 16),
            strip.text = element_text(size = 11),
            plot.title = element_text(size = 18, hjust = 0.5))
    # Save plot
    plot_prefix <- "exp3_plot2_tree_depth_Recent_freey_s"
    ggsave(p, filename = paste0(plot_directory, "/pdf_plots/", plot_prefix, s, ".pdf"), width = 10, height = 12.5, units = "in")
    ggsave(p, filename = paste0(plot_directory, "/png_plots/", plot_prefix, s, ".png"), width = 10, height = 12.5, units = "in")


    ## Plot 3: Ancient events. Smooth lines showing average values for each test statistic as the number of trees increases, faceted by tree number of taxa ##
    # Set dataset for plot
    plot_df <- s_df
    plot_df <- plot_df[plot_df$recombination_type == "Ancient", ]
    # Construct plot
    p <- ggplot(plot_df, aes(x = recombination_value, y = value, color = as.factor(tree_age))) +
      geom_smooth(method = "loess", alpha = 0.3, linewidth = 0, span = 0.75,
                  aes(ymin = ifelse(after_stat(ymin) < 0, 0, after_stat(ymin)), ymax = ifelse(after_stat(ymax) > 1, 1, after_stat(ymax)))) +
      stat_smooth(method = "loess", geom = "line", linewidth = 1.1, alpha = 1, span = 0.75) +
      facet_grid(var_label~num_taxa, scales = "fixed", labeller = label_parsed) +
      scale_x_continuous(name = "Proportion of recombinant DNA", breaks = seq(0,0.5, 0.25), labels = seq(0,0.5, 0.25), minor_breaks = seq(0,0.5, 0.05)) +
      scale_y_continuous(name = "Test statistic value", limits = c(0,1.10), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1), oob=scales::rescale_none) +
      scale_color_viridis_d(direction = -1, option = "C") +
      guides(color = guide_legend(title = "Tree depth\n(coalescent units)")) +
      labs(title = "Number of taxa") +
      theme_bw() +
      theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
            axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 14),
            legend.title = element_text(size = 18, hjust = 0.5), legend.text = element_text(size = 16),
            strip.text = element_text(size = 11),
            plot.title = element_text(size = 18, hjust = 0.5))
    # Save plot
    plot_prefix <- "exp3_plot3_num_taxa_Ancient_s"
    ggsave(p, filename = paste0(plot_directory, "/pdf_plots/", plot_prefix, s, ".pdf"), width = 10, height = 12.5, units = "in")
    ggsave(p, filename = paste0(plot_directory, "/png_plots/", plot_prefix, s, ".png"), width = 10, height = 12.5, units = "in")

    # Construct plot with points
    p <- ggplot(plot_df, aes(x = recombination_value, y = value, color = as.factor(tree_age))) +
      geom_point(alpha = 0.4) +
      facet_grid(var_label~num_taxa, scales = "fixed", labeller = label_parsed) +
      scale_x_continuous(name = "Proportion of recombinant DNA", breaks = seq(0,0.5, 0.25), labels = seq(0,0.5, 0.25), minor_breaks = seq(0,0.5, 0.05)) +
      scale_y_continuous(name = "Test statistic value", limits = c(0,1.10), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1), oob=scales::rescale_none) +
      scale_color_viridis_d(direction = -1, option = "C") +
      guides(color = guide_legend(title = "Tree depth\n(coalescent units)")) +
      labs(title = "Number of taxa") +
      theme_bw() +
      theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
            axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 14),
            legend.title = element_text(size = 18, hjust = 0.5), legend.text = element_text(size = 16),
            strip.text = element_text(size = 11),
            plot.title = element_text(size = 18, hjust = 0.5))
    # Save plot
    plot_prefix <- "exp3_plot3_num_taxa_Ancient_points_s"
    ggsave(p, filename = paste0(plot_directory, "/pdf_plots/", plot_prefix, s, ".pdf"), width = 10, height = 12.5, units = "in")
    ggsave(p, filename = paste0(plot_directory, "/png_plots/", plot_prefix, s, ".png"), width = 10, height = 12.5, units = "in")



    ## Plot 4: Recent events. Smooth lines showing average values for each test statistic as the number of trees increases, faceted by tree number of taxa ##
    # Set dataset for plot
    plot_df <- s_df
    plot_df <- plot_df[plot_df$recombination_type == "Recent", ]
    # Construct plot
    p <- ggplot(plot_df, aes(x = recombination_value, y = value, color = as.factor(tree_age))) +
      geom_smooth(method = "loess", alpha = 0.3, linewidth = 0, span = 0.75,
                  aes(ymin = ifelse(after_stat(ymin) < 0, 0, after_stat(ymin)), ymax = ifelse(after_stat(ymax) > 1, 1, after_stat(ymax)))) +
      stat_smooth(method = "loess", geom = "line", linewidth = 1.1, alpha = 1, span = 0.75) +
      facet_grid(var_label~num_taxa, scales = "fixed", labeller = label_parsed) +
      scale_x_continuous(name = "Proportion of recombinant DNA", breaks = seq(0,0.5, 0.1), labels = seq(0,0.5, 0.1), minor_breaks = seq(0,0.5, 0.05)) +
      scale_y_continuous(name = "Test statistic value", limits = c(0,1.10), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1), oob=scales::rescale_none) +
      scale_color_viridis_d(direction = -1, option = "C") +
      guides(color = guide_legend(title = "Tree depth\n(coalescent units)")) +
      labs(title = "Number of taxa") +
      theme_bw() +
      theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
            axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 14),
            legend.title = element_text(size = 18, hjust = 0.5), legend.text = element_text(size = 16),
            strip.text = element_text(size = 11),
            plot.title = element_text(size = 18, hjust = 0.5))
    # Save plot
    plot_prefix <- "exp3_plot4_num_taxa_Recent_s"
    ggsave(p, filename = paste0(plot_directory, "/pdf_plots/", plot_prefix, s, ".pdf"), width = 10, height = 12.5, units = "in")
    ggsave(p, filename = paste0(plot_directory, "/png_plots/", plot_prefix, s, ".png"), width = 10, height = 12.5, units = "in")

    # Construct plot with points
    p <- ggplot(plot_df, aes(x = recombination_value, y = value, color = as.factor(tree_age))) +
      geom_point(alpha = 0.4) +
      facet_grid(var_label~num_taxa, scales = "fixed", labeller = label_parsed) +
      scale_x_continuous(name = "Proportion of recombinant DNA", breaks = seq(0,0.5, 0.1), labels = seq(0,0.5, 0.1), minor_breaks = seq(0,0.5, 0.05)) +
      scale_y_continuous(name = "Test statistic value", limits = c(0,1.10), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1), oob=scales::rescale_none) +
      scale_color_viridis_d(direction = -1, option = "C") +
      guides(color = guide_legend(title = "Tree depth\n(coalescent units)")) +
      labs(title = "Number of taxa") +
      theme_bw() +
      theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
            axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 14),
            legend.title = element_text(size = 18, hjust = 0.5), legend.text = element_text(size = 16),
            strip.text = element_text(size = 11),
            plot.title = element_text(size = 18, hjust = 0.5))
    # Save plot
    plot_prefix <- "exp3_plot4_num_taxa_Recent_points_s"
    ggsave(p, filename = paste0(plot_directory, "/pdf_plots/", plot_prefix, s, ".pdf"), width = 10, height = 12.5, units = "in")
    ggsave(p, filename = paste0(plot_directory, "/png_plots/", plot_prefix, s, ".png"), width = 10, height = 12.5, units = "in")
  } # end iterating through substitution rates

} # end plotting experiment 3


