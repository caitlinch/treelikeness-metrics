# /caitlinch/treelikeness-metrics/03_data_analysis.R
# Caitlin Cherryh 2022

# This program takes results from applying various treelikeness tests and performs data analysis



#### 1. Set parameters ####
# data_directory          <- Directory where alignments will be saved/treelikeness metrics will be run.
# output_directory        <- Directory for output of data analysis
# repo_directory          <- Location of caitlinch/treelikeness-metrics github repository (for access to functions)

data_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/01_results/"
output_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/02_data_analysis/"
repo_directory <- "/Users/caitlincherryh/Documents/Repositories/treelikeness_metrics/"



#### 2. Prepare analyses ####
# Open packages
library(reshape2)
library(ggplot2)
library(patchwork)

# List all data files
data_files <- paste0(data_directory, list.files(data_directory))

# Create folder for plots
plot_directory <- paste0(output_directory, "plots/")
if (dir.exists(plot_directory) == FALSE){dir.create(plot_directory)}



#### 3. Prepare data from Experiment 1 for plotting ####
# Open data file from Experiment 1 as a dataframe
exp1_data_file <- grep("exp1", grep("treelikeness_metrics_collated_results", data_files, value = TRUE), value = TRUE)
exp1_df <- read.csv(file = exp1_data_file, stringsAsFactors = FALSE)
exp1_df$mean_TIGER_value <- as.numeric(exp1_df$mean_TIGER_value)

# Remove columns you don't want for plotting
exp1_wide_df <- exp1_df[, c("row_id", "uid", "num_taxa", "num_trees", "tree_depth", 
                            "tree_proportion", "Cunningham_test", "mean_delta_plot_value", 
                            "LM_proportion_resolved_quartets", "mean_Q_residual", 
                            "sCF_mean", "mean_TIGER_value")]

# Melt exp1_wide_df for better plotting
exp1_long_df <- melt(exp1_wide_df, id.vars = c("row_id", "uid", "num_taxa", "num_trees", "tree_depth"))

# Convert sCF values to decimal from percentage
exp1_wide_df$sCF_mean <- exp1_wide_df$sCF_mean / 100

# Transform the Network Treelikeness Test results into more plottable format
# Make a table of all possible parameter values for the network treelikeness test
ntlt_params <- expand.grid("num_taxa" = unique(exp1_df$num_taxa), "num_trees" = unique(exp1_df$num_trees), "tree_depth" = unique(exp1_df$tree_depth))
# Calculate proportion of treelike alignments for each set of parameter values
prop_tl_results <- unlist(lapply(1:nrow(ntlt_params), reformat.network.treelikeness.test.results, params_df = ntlt_params, results_df = exp1_df))
# Add columns to match the exp1_long_df
ntlt_params$row_id <- rep(NA, length(prop_tl_results))
ntlt_params$uid <- rep(NA, length(prop_tl_results))
ntlt_params$value <- prop_tl_results
ntlt_params$variable <- "NetworkTreelikenessTest"
# Restructure the dataframe to match the exp1_long_df
ntlt_params <- ntlt_params[,c(names(exp1_long_df))]
# Bind to the exp1_long_df
exp1_long_df <- rbind(exp1_long_df, ntlt_params)

# Add fancy labels for facets
exp1_long_df$var_label <- factor(exp1_long_df$variable, 
                                 levels = c("tree_proportion", "Cunningham_test", "mean_delta_plot_value", 
                                            "LM_proportion_resolved_quartets","NetworkTreelikenessTest",
                                            "mean_Q_residual", "sCF_mean", "mean_TIGER_value"), 
                                 ordered = TRUE, 
                                 labels = c(expression(atop("Tree","proportion")), expression(atop("Cunningham","metric")), 
                                            expression(paste('Mean ', delta["q"])), expression(atop("Proportion","resolved quartets")),
                                            expression(atop("Proportion","treelike alignments")), expression(atop("Mean", "Q-Residual value")), 
                                            expression(atop("Mean", "sCF value")), expression(atop("Mean","TIGER value"))) )

#### 4. Plot data from Experiment 1 ####
## Plot 1: Smooth lines showing average values for each test statistic as the number of trees increases, faceted by tree depth ##
# Set dataset for plot
plot_df <- exp1_long_df
# Set log10 minor breaks for x axis
x_axis_minor_breaks <-  unique(c(seq(1, 10, 1), seq(10, 100, 10), seq(100, 1000, 100), seq(1000, 10000, 1000)))
# Construct plot
p <- ggplot(plot_df, aes(x = num_trees, y = value, color = as.factor(num_taxa))) + 
  geom_smooth() + 
  facet_grid(var_label~tree_depth, scales = "fixed", labeller = label_parsed) +
  scale_x_log10( minor_breaks = x_axis_minor_breaks) +
  scale_y_continuous(name = "Test statistic value", limits = c(0,1.10), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_color_viridis_d(direction = -1) +
  guides(color = guide_legend(title = "Number of\ntaxa")) +
  labs(x = expression("Number of trees ("*log[10]*" scale)")) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 18), legend.text = element_text(size = 16),
        strip.text = element_text(size = 11))
# Save plot
plot_path <- paste0(plot_directory, "exp1_plot1_main.figure_tree_depth.pdf")
ggsave(p, filename = plot_path, width = 10, height = 12.5, units = "in")

## Plot 2: Smooth lines showing average values for each test statistic as the number of trees increases, faceted by tree number of taxa ##
# Set dataset for plot
plot_df <- exp1_long_df
# Set log10 minor breaks for x axis
x_axis_minor_breaks <-  unique(c(seq(1, 10, 1), seq(10, 100, 10), seq(100, 1000, 100), seq(1000, 10000, 1000)))
# Construct plot
p <- ggplot(plot_df, aes(x = num_trees, y = value, color = as.factor(tree_depth))) + 
  geom_smooth() + 
  facet_grid(var_label~num_taxa, scales = "fixed", labeller = label_parsed) +
  scale_x_log10( minor_breaks = x_axis_minor_breaks) +
  scale_y_continuous(name = "Test statistic value", limits = c(0,1.10), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_color_viridis_d(direction = -1) +
  guides(color = guide_legend(title = "Tree depth\n(substitutions\nper site)")) +
  labs(x = expression("Number of trees ("*log[10]*" scale)")) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 18), legend.text = element_text(size = 16),
        strip.text = element_text(size = 11))
# Save plot
plot_path <- paste0(plot_directory, "exp1_plot2_main.figure_num_taxa.pdf")
ggsave(p, filename = plot_path, width = 10, height = 12, units = "in")






#### 4. Plot data from Experiment 1 ####
## Plot 1: Smooth lines showing average values for each test statistic as the number of trees increases, colored by number of taxa ##
# Set dataset for plot. Extract only tree depth 0.1
plot_df <- exp1_long_df[exp1_long_df$tree_depth == 0.1, ]
# Construct plot
p <- ggplot(plot_df, aes(x = num_trees, y = value, color = as.factor(num_taxa))) + 
  geom_smooth() + 
  facet_wrap(~var_label, scales = "fixed", labeller = label_parsed, nrow = 3 ) +
  scale_x_log10(name = "Number of trees present in the alignment") +
  scale_y_continuous(name = "Test statistic value", limits = c(0,1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_color_viridis_d(direction = -1) +
  guides(color = guide_legend(title = "Number of\ntaxa")) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 18), legend.text = element_text(size = 14),
        strip.text = element_text(size = 13))
# Save plot
plot_path <- paste0(plot_directory, "exp1_plot1_main_figure_tree-depth-0.1.pdf")
ggsave(p, filename = plot_path)

## Plot 2: Smooth lines showing average values for each test statistic as the number of trees increases, faceted by tree depth ##
# Set dataset for plot
plot_df <- exp1_long_df
# Construct plot
p <- ggplot(plot_df, aes(x = num_trees, y = value, color = as.factor(num_taxa))) + 
  geom_smooth() + 
  facet_grid(var_label~tree_depth, scales = "fixed", labeller = label_parsed) +
  scale_x_log10(name = "Number of trees present in the alignment") +
  scale_y_continuous(name = "Test statistic value", limits = c(0,1), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_color_viridis_d(direction = -1) +
  guides(color = guide_legend(title = "Number of\ntaxa")) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 18), legend.text = element_text(size = 14),
        strip.text = element_text(size = 13))
# Save plot
plot_path <- paste0(plot_directory, "exp1_plot2_supplementary_tree_depth.pdf")
ggsave(p, filename = plot_path, height = 14, width = 11, units = "in")

## Plot 3: Smooth lines showing average values for each test statistic as the number of trees increases, faceted by number of taxa ##
# Set dataset for plot
plot_df <- exp1_long_df
# Construct plot
p <- ggplot(plot_df, aes(x = num_trees, y = value, color = as.factor(tree_depth))) + 
  geom_smooth() + 
  facet_grid(var_label~num_taxa, scales = "fixed", labeller = label_parsed) +
  scale_x_log10(name = "Number of trees present in the alignment") +
  scale_y_continuous(name = "Test statistic value", limits = c(0,1), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_color_viridis_d(direction = -1, option = "C") +
  guides(color = guide_legend(title = "Tree depth\n(substitutions\nper site)")) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 16), legend.text = element_text(size = 14),
        strip.text = element_text(size = 13))
# Save plot
plot_path <- paste0(plot_directory, "exp1_plot3_supplementary_num_taxa.pdf")
ggsave(p, filename = plot_path, height = 14, width = 11, units = "in")

## Plot 4: Bar chart for Network Treelikeness test as the number of trees increases ##
# Set dataset for plot
plot_df <- exp1_df
# Construct plot
p <- ggplot(plot_df, aes(x = NetworkTreelikenessTest)) +
  geom_bar() +
  facet_wrap(~num_trees) +
  scale_y_continuous(name = "Count", breaks = seq(0,200,50), labels = seq(0,200,50), 
                     minor_breaks = seq(0,200,25), limits = c(0, 150)) +
  scale_x_discrete(name = "Network treelikeness test result") +
  labs(title = "Number of taxa") + 
  theme_bw() + 
  theme(axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 13, angle = 30, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 13),
        strip.text = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5))
# Save plot
plot_path <- paste0(plot_directory, "exp1_plot4_main_figure_NetworkTreelikenessTest.pdf")
ggsave(p, filename = plot_path)

## Plot 3: Smooth lines showing average values for each test statistic as the number of trees increases, faceted by number of taxa ##
# Set dataset for plot
plot_df <- exp1_long_df
# Construct plot
p <- ggplot(plot_df, aes(x = num_trees, y = value, color = as.factor(tree_depth))) + 
  geom_smooth() + 
  facet_grid(var_label~num_taxa, scales = "fixed", labeller = label_parsed) +
  scale_x_log10(name = "Number of trees present in the alignment") +
  scale_y_continuous(name = "Test statistic value", limits = c(0,1), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_color_viridis_d(direction = -1, option = "C") +
  guides(color = guide_legend(title = "Tree depth\n(substitutions\nper site)")) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 16), legend.text = element_text(size = 14),
        strip.text = element_text(size = 13))
# Save plot
plot_path <- paste0(plot_directory, "exp1_plot3_supplementary_num_taxa.pdf")
ggsave(p, filename = plot_path, height = 14, width = 11, units = "in")

## Plot 4: Bar chart for Network Treelikeness test as the number of trees increases, faceted by tree depth ##
# Set dataset for plot
plot_df <- exp1_df
# Construct plot
p <- ggplot(plot_df, aes(x = NetworkTreelikenessTest)) +
  geom_bar() +
  facet_wrap(~num_trees) +
  scale_y_continuous(name = "Count", breaks = seq(0,200,50), labels = seq(0,200,50), 
                     minor_breaks = seq(0,200,25), limits = c(0, 150)) +
  scale_x_discrete(name = "Network treelikeness test result") +
  labs(title = "Number of taxa") + 
  theme_bw() + 
  theme(axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 13, angle = 30, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 13),
        strip.text = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5))
# Save plot
plot_path <- paste0(plot_directory, "exp1_plot4_main_figure_NetworkTreelikenessTest.pdf")
ggsave(p, filename = plot_path)

## Plot 5: Bar chart for Network Treelikeness test as the number of trees increases, faceted by tree depth ##
# Set dataset for plot
plot_df <- exp1_df
# Construct plot
p <- ggplot(plot_df, aes(x = NetworkTreelikenessTest)) +
  geom_bar() +
  facet_grid(num_trees ~ tree_depth) +
  scale_y_continuous(name = "Count", breaks = seq(0,60,30), labels = seq(0,60,30), 
                     minor_breaks = seq(0,60,10), limits = c(0, 60)) +
  scale_x_discrete(name = "Network treelikeness test result") +
  labs(title = "Tree depth") + 
  theme_bw() + 
  theme(axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 11, angle = 30, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 11),
        strip.text = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5))
# Save plot
plot_path <- paste0(plot_directory, "exp1_plot5_supplementary_NetworkTreelikenessTest_tree-depth.pdf")
ggsave(p, filename = plot_path, height = 20, width = 11, units = "in")

## Plot 6: Bar chart for Network Treelikeness test as the number of trees increases, faceted by number of taxa ##
# Set dataset for plot
plot_df <- exp1_df
# Construct plot
p <- ggplot(plot_df, aes(x = NetworkTreelikenessTest)) +
  geom_bar() +
  facet_grid(num_trees ~ num_taxa) +
  scale_y_continuous(name = "Count", breaks = seq(0,40,20), labels = seq(0,40,20), 
                     minor_breaks = seq(0,40,10), limits = c(0, 40)) +
  scale_x_discrete(name = "Network treelikeness test result") +
  labs(title = "Number of taxa") + 
  theme_bw() + 
  theme(axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 11, angle = 30, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 11),
        strip.text = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5))
# Save plot
plot_path <- paste0(plot_directory, "exp1_plot6_supplementary_NetworkTreelikenessTest_num-taxa.pdf")
ggsave(p, filename = plot_path, height = 20, width = 11, units = "in")


