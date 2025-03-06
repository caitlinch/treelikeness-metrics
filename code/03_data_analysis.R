# caitlinch/treelikeness-metrics/code/03_data_analysis.R
# Caitlin Cherryh 2025

# This program takes results from applying various treelikeness tests and
# performs data analysis/plotting

# Use the treelikeness-metrics RProject for this script



#### 1. Set parameters ####
# data_directory <- Directory containing result csv files
# plot_directory <- Directory for output of data analysis/plots
# save.png <- TRUE to output png plots (1-2 MB per figure)
data_directory <- "output/"
plot_directory <- "plot/"
save.png <- TRUE



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
if (save.png){
  ggsave(exp1_plot1,
         filename = paste0(plot_directory, "exp1_plot1_median_ribbon.png"),
         width = 10,
         height = 13,
         units = "in"
  )
}


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
if (save.png){
  ggsave(exp1_plot2,
         filename = paste0(plot_directory, "exp1_plot2_median_ribbon_freeY.png"),
         width = 12,
         height = 14,
         units = "in"
  )
}

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
if (save.png){
  ggsave(exp1_plot3,
         filename = paste0(plot_directory, "exp1_plot3_points.png"),
         width = 10,
         height = 13,
         units = "in"
  )
}



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
exp3_df$sCF_mean <- exp3_df$sCF_mean / 100

# Convert mean tiger value to numeric
# Note: 4% (150/3750) of simulations for experiment 1 do not have a TIGER value (TIGER failed to run in these cases)
#       Therefore converting to numeric will coerce these values (i.e. mean_TIGER_value = "no_TIGER_run") to NA
exp3_df$mean_TIGER_value <- as.numeric(exp3_df$mean_TIGER_value)

# Remove columns you don't want for plotting
nonbinary_metrics <- c(
  "tree_proportion",
  "Cunningham_test",
  "mean_delta_plot_value",
  "LM_proportion_resolved_quartets",
  "mean_Q_residual",
  "sCF_mean",
  "mean_TIGER_value"
)
exp3_df <- exp3_df[, c(
  "uid",
  "row_id",
  "num_reps",
  "num_taxa",
  "num_trees",
  "tree_age",
  "tree_depth_subspersite",
  "recombination_value",
  "recombination_type",
  "speciation_rate",
  nonbinary_metrics,
  "NetworkTreelikenessTest"
)]

# Create data frame of summary statistics for the non-binary metrics
exp3_summary_df <- as.data.frame(do.call(
  rbind,
  lapply(nonbinary_metrics, exp3.metric.statistics.wrapper, df = exp3_df)
))

# Reformat and bind the Network Treelikeness Test results
exp3_summary_df <- rbind(exp3_summary_df,
                         exp3.process.NetworkTreelikenessTest(exp3_df))

# Format the experiment 3 results for pretty plots
exp3_df <- exp3.format.dataframe(exp3_df)
exp3_summary_df <- exp3.format.summary.dataframe(exp3_summary_df)



#### 6. Plot data from Experiment 3 ####
# Create list of plots to generate for the different experimental conditions
plot_list <-
  list(
    c("recombination_type" = "Ancient", "speciation_rate" = 0.1),
    c("recombination_type" = "Ancient", "speciation_rate" = 1.0),
    c("recombination_type" = "Recent", "speciation_rate" = 0.1),
    c("recombination_type" = "Recent", "speciation_rate" = 1.0)
  )

## Plot 1:
# Lines showing median values for each test statistic as the number of trees
# increases, faceted by tree depth
#   - Log10 x-axis
#   - Free y-axis for each test statistic
for (i in 1:length(plot_list)){
  plot_properties <- plot_list[[i]]
  print(paste0(
    "Experiment 3 - Plot 1 - ",
    plot_properties[["recombination_type"]],
    " - s",
    plot_properties[["speciation_rate"]]
  ))
  if (plot_properties[["recombination_type"]] == "Ancient"){
    plot_palette = scales::viridis_pal(direction = -1)(5)[2:5]
  } else if (plot_properties[["recombination_type"]] == "Recent"){
    plot_palette = scales::viridis_pal(direction = -1)(5)
  }
  plot_title = paste0(plot_properties[["recombination_type"]],
                      " introgression; speciation rate = ",
                      plot_properties[["speciation_rate"]])

  plot_df <- exp3_summary_df[which(
    exp3_summary_df$recombination_type == plot_properties[["recombination_type"]] &
      exp3_summary_df$speciation_rate == as.numeric(plot_properties[["speciation_rate"]])
  ), ]
  exp3_plot1 <- ggplot(plot_df) +
    geom_line(aes(
      x = recombination_value,
      y = median,
      color = as.factor(num_taxa)
    )) +
    geom_ribbon(aes(
      x = recombination_value,
      ymin = minimum,
      ymax = maximum,
      fill = as.factor(num_taxa)
    ),
    alpha = 0.15) +
    facet_grid(var_label ~ tree_age, scales = "free", labeller = label_parsed) +
    scale_x_continuous(
      name = "Proportion of recombinant DNA",
      breaks = seq(0, 0.5, 0.1),
      labels = seq(0, 0.5, 0.1),
      minor_breaks = seq(0, 0.5, 0.05)
    ) +
    scale_y_continuous(name = "Median test statistic value") +
    scale_color_manual(values = plot_palette) +
    scale_fill_manual(values = plot_palette) +
    guides(
      color = guide_legend(title = "Number of\ntaxa"),
      fill = guide_legend(title = "Number of\ntaxa")
    ) +
    labs(title = plot_title, subtitle = "Tree age (Million years)") +
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
      legend.title = element_text(size = 18, hjust = 0.5),
      legend.text = element_text(size = 16),
      strip.text = element_text(size = 11),
      plot.title = element_text(
        hjust = 0.5,
        size = 25,
        margin = margin(t = 0, r = 0, b = 15, l = 0)
      ),
      plot.subtitle = element_text(size = 18, hjust = 0.5)
    )

  plot_file_suffix <- paste0(
    plot_properties[["recombination_type"]],
    "_s",
    plot_properties[["speciation_rate"]]
  )
  ggsave(exp3_plot1,
         filename = paste0(plot_directory, "exp3_plot1_", plot_file_suffix, "_median_ribbon.pdf"),
         width = 10,
         height = 13,
         units = "in"
  )
  if (save.png){
    ggsave(exp3_plot1,
           filename = paste0(plot_directory, "exp3_plot1_", plot_file_suffix, "_median_ribbon.png"),
           width = 10,
           height = 13,
           units = "in"
    )
  }
}

## Plot 2:
# Lines showing median values for each test statistic as the number of trees
# increases, faceted by tree depth
#   - Log10 x-axis
#   - Free y-axis for each test statistic and each tree depth
for (i in 1:length(plot_list)){
  plot_properties <- plot_list[[i]]
  print(paste0(
    "Experiment 3 - Plot 2 - ",
    plot_properties[["recombination_type"]],
    " - s",
    plot_properties[["speciation_rate"]]
  ))
  if (plot_properties[["recombination_type"]] == "Ancient"){
    plot_palette = scales::viridis_pal(direction = -1)(5)[2:5]
  } else if (plot_properties[["recombination_type"]] == "Recent"){
    plot_palette = scales::viridis_pal(direction = -1)(5)
  }
  plot_title = paste0(plot_properties[["recombination_type"]],
                      " introgression; speciation rate = ",
                      plot_properties[["speciation_rate"]])
  plot_df <- exp3_summary_df[which(
    exp3_summary_df$recombination_type == plot_properties[["recombination_type"]] &
      exp3_summary_df$speciation_rate == as.numeric(plot_properties[["speciation_rate"]])
  ), ]

  exp3_plot2_p1 <-
    ggplot(plot_df[which(plot_df$tree_age == 5), ]) +
    geom_line(aes(
      x = recombination_value,
      y = median,
      color = as.factor(num_taxa)
    )) +
    geom_ribbon(aes(
      x = recombination_value,
      ymin = minimum,
      ymax = maximum,
      fill = as.factor(num_taxa)
    ),
    alpha = 0.15) +
    facet_grid(var_label ~ tree_age, scales = "free", labeller = label_parsed) +
    scale_x_continuous(
      name = "Proportion of recombinant DNA",
      breaks = seq(0, 0.5, 0.1),
      labels = seq(0, 0.5, 0.1),
      minor_breaks = seq(0, 0.5, 0.05)
    ) +
    scale_y_continuous(name = "Median test statistic value") +
    scale_color_manual(values = plot_palette) +
    scale_fill_manual(values = plot_palette) +
    guides(
      color = guide_legend(title = "Number of\ntaxa"),
      fill = guide_legend(title = "Number of\ntaxa")
    ) +
    labs(title = plot_title, subtitle = "Tree age (Million years)") +
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

  exp3_plot2_p2 <-
    ggplot(plot_df[which(plot_df$tree_age == 50), ]) +
    geom_line(aes(
      x = recombination_value,
      y = median,
      color = as.factor(num_taxa)
    )) +
    geom_ribbon(aes(
      x = recombination_value,
      ymin = minimum,
      ymax = maximum,
      fill = as.factor(num_taxa)
    ),
    alpha = 0.15) +
    facet_grid(var_label ~ tree_age, scales = "free", labeller = label_parsed) +
    scale_x_continuous(
      name = "Proportion of recombinant DNA",
      breaks = seq(0, 0.5, 0.1),
      labels = seq(0, 0.5, 0.1),
      minor_breaks = seq(0, 0.5, 0.05)
    ) +
    scale_y_continuous(name = "Median test statistic value") +
    scale_color_manual(values = plot_palette) +
    scale_fill_manual(values = plot_palette) +
    guides(
      color = guide_legend(title = "Number of\ntaxa"),
      fill = guide_legend(title = "Number of\ntaxa")
    ) +
    labs(title = plot_title, subtitle = "Tree age (Million years)") +
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

  exp3_plot2_p3 <-
    ggplot(plot_df[which(plot_df$tree_age == 500), ]) +
    geom_line(aes(
      x = recombination_value,
      y = median,
      color = as.factor(num_taxa)
    )) +
    geom_ribbon(aes(
      x = recombination_value,
      ymin = minimum,
      ymax = maximum,
      fill = as.factor(num_taxa)
    ),
    alpha = 0.15) +
    facet_grid(var_label ~ tree_age, scales = "free", labeller = label_parsed) +
    scale_x_continuous(
      name = "Proportion of recombinant DNA",
      breaks = seq(0, 0.5, 0.1),
      labels = seq(0, 0.5, 0.1),
      minor_breaks = seq(0, 0.5, 0.05)
    ) +
    scale_y_continuous(name = "Median test statistic value") +
    scale_color_manual(values = plot_palette) +
    scale_fill_manual(values = plot_palette) +
    guides(
      color = guide_legend(title = "Number of\ntaxa"),
      fill = guide_legend(title = "Number of\ntaxa")
    ) +
    labs(title = plot_title, subtitle = "Tree age (Million years)") +
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

  exp3_plot2 <-
    exp3_plot2_p1 + exp3_plot2_p2 + exp3_plot2_p3

  plot_file_suffix <- paste0(
    plot_properties[["recombination_type"]],
    "_s",
    plot_properties[["speciation_rate"]]
  )
  ggsave(exp3_plot2,
         filename = paste0(plot_directory, "exp3_plot2_", plot_file_suffix, "_median_ribbon_freeY.pdf"),
         width = 10,
         height = 13,
         units = "in"
  )
  if (save.png){
    ggsave(exp3_plot2,
           filename = paste0(plot_directory, "exp3_plot2_", plot_file_suffix, "_median_ribbon_freeY.png"),
           width = 10,
           height = 13,
           units = "in"
    )
  }
}

## Plot 3:
# Dot plot of raw data
#   - Log10 x-axis
#   - Free y-axis for each column
# Note: 150 NAs from TIGER, 2 values < 0 from Cunningham metric
for (i in 1:length(plot_list)){
  plot_properties <- plot_list[[i]]
  print(paste0(
    "Experiment 3 - Plot 3 - ",
    plot_properties[["recombination_type"]],
    " - s",
    plot_properties[["speciation_rate"]]
  ))
  if (plot_properties[["recombination_type"]] == "Ancient"){
    plot_palette = scales::viridis_pal(direction = -1)(5)[2:5]
  } else if (plot_properties[["recombination_type"]] == "Recent"){
    plot_palette = scales::viridis_pal(direction = -1)(5)
  }
  plot_title = paste0(plot_properties[["recombination_type"]],
                      " introgression; speciation rate = ",
                      plot_properties[["speciation_rate"]])
  plot_df <- exp3_df[which(
    exp3_df$recombination_type == plot_properties[["recombination_type"]] &
      exp3_df$speciation_rate == as.numeric(plot_properties[["speciation_rate"]])
  ), ]

  exp3_plot3 <- ggplot(plot_df) +
    geom_point(aes(
      x = recombination_value,
      y = value,
      color = as.factor(num_taxa)
    )) +
    facet_grid(var_label ~ tree_age, scales = "free", labeller = label_parsed) +
    scale_x_continuous(
      name = "Proportion of recombinant DNA",
      breaks = seq(0, 0.5, 0.1),
      labels = seq(0, 0.5, 0.1),
      minor_breaks = seq(0, 0.5, 0.05)
    ) +
    scale_y_continuous(name = "Median test statistic value") +
    scale_color_manual(values = plot_palette) +
    scale_fill_manual(values = plot_palette) +
    guides(
      color = guide_legend(title = "Number of\ntaxa"),
      fill = guide_legend(title = "Number of\ntaxa")
    ) +
    labs(title = plot_title, subtitle = "Tree age (Million years)") +
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
      legend.title = element_text(size = 18, hjust = 0.5),
      legend.text = element_text(size = 16),
      strip.text = element_text(size = 11),
      plot.title = element_text(
        hjust = 0.5,
        size = 25,
        margin = margin(t = 0, r = 0, b = 15, l = 0)
      ),
      plot.subtitle = element_text(size = 18, hjust = 0.5)
    )

  plot_file_suffix <- paste0(
    plot_properties[["recombination_type"]],
    "_s",
    plot_properties[["speciation_rate"]]
  )
  ggsave(exp3_plot3,
         filename = paste0(plot_directory, "exp3_plot3_", plot_file_suffix, "_points.pdf"),
         width = 10,
         height = 13,
         units = "in"
  )
  if (save.png){
    ggsave(exp3_plot3,
           filename = paste0(plot_directory, "exp3_plot3_", plot_file_suffix, "_points.png"),
           width = 10,
           height = 13,
           units = "in"
    )
  }
}


