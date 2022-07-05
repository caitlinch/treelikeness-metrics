# /caitlinch/treelikeness-metrics/03_data_analysis.R
# Caitlin Cherryh 2022

# This program takes results from applying various treelikeness tests and performs data analysis



#### 1. Set parameters ####
# data_directory          <- Directory where alignments will be saved/treelikeness metrics will be run.
# output_directory        <- Directory for output of data analysis
# repo_directory          <- Location of caitlinch/treelikeness-metrics github repository (for access to functions).

data_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/01_results/"
output_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/02_data_analysis/"
repo_directory <- "/Users/caitlincherryh/Documents/Repositories/treelikeness_metrics/"


#### 2. Prepare analyses ####
# Open packages
library(dplyr) # for mutate
library(reshape2)
library(ggplot2)
library(patchwork)

# List all data files
data_files <- paste0(data_directory, list.files(data_directory))


# Create folder for plots
plot_directory <- paste0(output_directory, "plots/")
if (dir.exists(plot_directory) == FALSE){dir.create(plot_directory)}

#### 3. Analyse output from Experiment 1 ####
# Open data file from Experiment 1 as a dataframe
exp1_data_file <- grep("exp1", grep("treelikeness_metrics_results", data_files, value = TRUE), value = TRUE)
exp1_df <- read.csv(file = exp1_data_file, stringsAsFactors = FALSE)

## Compare ids with those in parameters csv for experiment 1 to determine if there are any incomplete/missing alignments











# Remove columns you don't want for plotting
exp1_wide_df <- exp1_df[, c("row_id", "uid", "num_taxa", "num_trees", "tree_depth", 
                            "tree_proportion", "Cunningham_test", "mean_delta_plot_value", 
                            "LM_proportion_resolved_quartets", "mean_Q_residual", 
                            "sCF_mean", "mean_TIGER_value")]
exp1_ntlt_df <-exp1_df[, c("row_id", "uid", "num_taxa", "num_trees", "tree_depth", 
                           "NetworkTreelikenessTest")]
# Melt exp1_wide_df for better plotting
exp1_long_df <- melt(exp1_wide_df, id.vars = c("row_id", "uid", "num_taxa", "num_trees", "tree_depth"))
# Add fancy labels for facets
exp1_long_df$var_label <- factor(exp1_long_df$variable, 
                                 levels = c("tree_proportion", "Cunningham_test", "mean_delta_plot_value", 
                                            "LM_proportion_resolved_quartets", "mean_Q_residual", "sCF_mean",
                                            "mean_TIGER_value"), 
                                 ordered = TRUE, 
                                 labels = c(expression("`Tree proportion`"), expression("`Cunningham metric`"), 
                                            expression(paste('Mean ', delta["q"])), expression("`Proportion resolved quartets`"),
                                            expression("`Mean Q-Residual value`"), expression("`Mean sCF value`"),
                                            expression("`Mean TIGER value`")) )


# Extract only tree depth 0.1 and taxa = 50
exp1_long_df <- exp1_long_df[exp1_long_df$tree_depth == 0.1 & exp1_long_df$num_taxa == 50, ]

# Plot 1: Faceted scatterplot
ggplot(exp1_long_df, aes(x = num_trees, y = value, color = factor(num_taxa))) + 
  geom_smooth() + 
  facet_wrap(~var_label, nrow = 2, scales = "free_y", labeller = label_parsed) +
  theme_bw()

ggplot(exp1_long_df, aes(x = num_trees, y = value)) + 
  geom_smooth() + 
  facet_grid(var_label ~ num_taxa, scales = "free_y", labeller = label_parsed) +
  theme_bw()


plot_path <- paste0(plot_directory, "exp1_Plot_1_faceted_scatterplot")
