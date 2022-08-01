# /caitlinch/treelikeness-metrics/code/04_data_analysis.R
# Caitlin Cherryh 2022

# This program takes results from applying various treelikeness tests and performs data analysis



#### 1. Set parameters ####
# data_directory          <- Directory where alignments will be saved/treelikeness metrics will be run.
# output_directory        <- Directory for output of data analysis
# repo_directory          <- Location of caitlinch/treelikeness-metrics github repository (for access to functions)

# plot_exp1 <- either TRUE to plot experiment 1 results, or FALSE to skip
# plot_exp2 <- either TRUE to plot experiment 2 results, or FALSE to skip
# plot_empirical <- either TRUE to plot empirical dataset results, or FALSE to skip

data_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/01_results/"
output_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/02_data_analysis/"
repo_directory <- "/Users/caitlincherryh/Documents/Repositories/treelikeness-metrics/"

plot_exp1 = FALSE
plot_exp2 = FALSE
plot_empirical = TRUE


#### 2. Prepare analyses ####
# Source function files
source(paste0(repo_directory, "code/func_data_analysis.R"))

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
if (plot_exp1 == TRUE){
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
  prop_tl_results <- unlist(lapply(1:nrow(ntlt_params), reformat.network.treelikeness.test.results.exp1, params_df = ntlt_params, results_df = exp1_df))
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
}



#### 4. Plot data from Experiment 1 ####
if (plot_exp1 == TRUE){
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
    scale_color_viridis_d(direction = -1, option = "C") +
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
}



#### 5. Prepare data from Experiment 2 for plotting ####
if (plot_exp2 == TRUE){
  # Open data file from Experiment 2 as a dataframe
  exp2_data_file <- grep("exp2", grep("treelikeness_metrics_collated_results", data_files, value = TRUE), value = TRUE)
  exp2_df <- read.csv(file = exp2_data_file, stringsAsFactors = FALSE)
  exp2_df$mean_TIGER_value <- as.numeric(exp2_df$mean_TIGER_value)
  
  # Convert sCF values to decimal from percentage
  exp2_df$sCF_mean <- exp2_df$sCF_mean / 100
  
  # # Remove columns you don't want for plotting
  # exp2_wide_df <- exp2_df[, c("row_id", "uid", "num_taxa", "num_trees", "tree_depth", 
  #                             "recombination_value", "recombination_type",
  #                             "tree_proportion", "Cunningham_test", "mean_delta_plot_value", 
  #                             "LM_proportion_resolved_quartets", "mean_Q_residual", 
  #                             "sCF_mean", "mean_TIGER_value")]
  # Remove columns you don't want for plotting - do not plot TIGER (fast TIGER was not run for exp2, too time consuming)
  exp2_wide_df <- exp2_df[, c("row_id", "uid", "num_taxa", "num_trees", "tree_depth", 
                              "recombination_value", "recombination_type",
                              "tree_proportion", "Cunningham_test", "mean_delta_plot_value", 
                              "LM_proportion_resolved_quartets", "mean_Q_residual", 
                              "sCF_mean")]
  
  # Melt exp2_wide_df for better plotting
  exp2_long_df <- melt(exp2_wide_df, id.vars = c("row_id", "uid", "num_taxa", "num_trees", "tree_depth", "recombination_value", "recombination_type"))
  
  
  # Transform the Network Treelikeness Test results into more plottable format
  # Make a table of all possible parameter values for the network treelikeness test
  ntlt_params <- expand.grid("num_taxa" = sort(unique(exp2_df$num_taxa)), "tree_depth" = sort(unique(exp2_df$tree_depth)), 
                             "recombination_value" = sort(unique(exp2_df$recombination_value)), "recombination_type" = unique(exp2_df$recombination_type))
  # Calculate proportion of treelike alignments for each set of parameter values
  prop_tl_results <- unlist(lapply(1:nrow(ntlt_params), reformat.network.treelikeness.test.results.exp2, params_df = ntlt_params, results_df = exp2_df))
  # Add columns to match the exp2_long_df
  ntlt_params$row_id <- rep(NA, length(prop_tl_results))
  ntlt_params$uid <- rep(NA, length(prop_tl_results))
  ntlt_params$num_trees <- 1000
  ntlt_params$value <- prop_tl_results
  ntlt_params$variable <- "NetworkTreelikenessTest"
  # Restructure the dataframe to match the exp2_long_df
  ntlt_params <- ntlt_params[,c(names(exp2_long_df)),]
  # Bind to the exp2_long_df
  exp2_long_df <- rbind(exp2_long_df, ntlt_params)
  
  # # Add fancy labels for facets
  # exp2_long_df$var_label <- factor(exp2_long_df$variable, 
  #                                  levels = c("tree_proportion", "Cunningham_test", "mean_delta_plot_value", 
  #                                             "LM_proportion_resolved_quartets","NetworkTreelikenessTest",
  #                                             "mean_Q_residual", "sCF_mean", "mean_TIGER_value"), 
  #                                  ordered = TRUE, 
  #                                  labels = c(expression(atop("Tree","proportion")), expression(atop("Cunningham","metric")), 
  #                                             expression(paste('Mean ', delta["q"])), expression(atop("Proportion","resolved quartets")),
  #                                             expression(atop("Proportion","treelike alignments")), expression(atop("Mean", "Q-Residual value")), 
  #                                             expression(atop("Mean", "sCF value")), expression(atop("Mean","TIGER value"))) )
  # Add fancy labels for facets -  - do not plot TIGER (fast TIGER was not run for exp2, too time consuming)
  exp2_long_df$var_label <- factor(exp2_long_df$variable, 
                                   levels = c("tree_proportion", "Cunningham_test", "mean_delta_plot_value", 
                                              "LM_proportion_resolved_quartets","NetworkTreelikenessTest",
                                              "mean_Q_residual", "sCF_mean"), 
                                   ordered = TRUE, 
                                   labels = c(expression(atop("Tree","proportion")), expression(atop("Cunningham","metric")), 
                                              expression(paste('Mean ', delta["q"])), expression(atop("Proportion","resolved quartets")),
                                              expression(atop("Proportion","treelike alignments")), expression(atop("Mean", "Q-Residual value")), 
                                              expression(atop("Mean", "sCF value")) ) )
}



#### 6. Plot data from Experiment 2 ####
if (plot_exp2 == TRUE){
  ## Plot 1: Ancient events. Smooth lines showing average values for each test statistic as the amount of recombination increases, faceted by tree depth ##
  # Set dataset for plot
  plot_df <- exp2_long_df
  plot_df <- plot_df[plot_df$recombination_type == "Ancient", ]
  # Construct plot
  p <- ggplot(plot_df, aes(x = recombination_value, y = value, color = as.factor(num_taxa))) + 
    geom_smooth() + 
    facet_grid(var_label~tree_depth, scales = "fixed", labeller = label_parsed) +
    scale_x_continuous(name = "Proportion of recombinant DNA", breaks = seq(0,0.5, 0.1), labels = seq(0,0.5, 0.1), minor_breaks = seq(0,0.5, 0.05)) +
    scale_y_continuous(name = "Test statistic value", limits = c(0,1.10), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_color_viridis_d(direction = -1) +
    guides(color = guide_legend(title = "Number of\ntaxa")) +
    labs(title = "Tree depth (coalescent units)") +
    theme_bw() +
    theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
          axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 14),
          legend.title = element_text(size = 18, hjust = 0.5), legend.text = element_text(size = 16),
          strip.text = element_text(size = 11),
          plot.title = element_text(size = 18, hjust = 0.5))
  # Save plot
  plot_path <- paste0(plot_directory, "exp2_plot1_main.figure_tree_depth_Ancient.pdf")
  ggsave(p, filename = plot_path, width = 10, height = 12.5, units = "in")
  
  ## Plot 2: Recent events. Smooth lines showing average values for each test statistic as the amount of recombination increases, faceted by tree depth ##
  # Set dataset for plot
  plot_df <- exp2_long_df
  plot_df <- plot_df[plot_df$recombination_type == "Recent", ]
  # Construct plot
  p <- ggplot(plot_df, aes(x = recombination_value, y = value, color = as.factor(num_taxa))) + 
    geom_smooth() + 
    facet_grid(var_label~tree_depth, scales = "fixed", labeller = label_parsed) +
    scale_x_continuous(name = "Proportion of recombinant DNA", breaks = seq(0,0.5, 0.1), labels = seq(0,0.5, 0.1), minor_breaks = seq(0,0.5, 0.05)) +
    scale_y_continuous(name = "Test statistic value", limits = c(0,1.10), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_color_viridis_d(direction = -1) +
    guides(color = guide_legend(title = "Number of\ntaxa")) +
    labs(title = "Tree depth (coalescent units)") +
    theme_bw() +
    theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
          axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 14),
          legend.title = element_text(size = 18, hjust = 0.5), legend.text = element_text(size = 16),
          strip.text = element_text(size = 11),
          plot.title = element_text(size = 18, hjust = 0.5))
  # Save plot
  plot_path <- paste0(plot_directory, "exp2_plot2_main.figure_tree_depth_Recent.pdf")
  ggsave(p, filename = plot_path, width = 10, height = 12.5, units = "in")
  
  ## Plot 3: Ancient events. Smooth lines showing average values for each test statistic as the number of trees increases, faceted by tree number of taxa ##
  # Set dataset for plot
  plot_df <- exp2_long_df
  plot_df <- plot_df[plot_df$recombination_type == "Ancient", ]
  # Construct plot
  p <- ggplot(plot_df, aes(x = recombination_value, y = value, color = as.factor(tree_depth))) + 
    geom_smooth() + 
    facet_grid(var_label~num_taxa, scales = "fixed", labeller = label_parsed) +
    scale_x_continuous(name = "Proportion of recombinant DNA", breaks = seq(0,0.5, 0.25), labels = seq(0,0.5, 0.25), minor_breaks = seq(0,0.5, 0.05)) +
    scale_y_continuous(name = "Test statistic value", limits = c(0,1.10), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1)) +
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
  plot_path <- paste0(plot_directory, "exp2_plot3_main.figure_num_taxa_Ancient.pdf")
  ggsave(p, filename = plot_path, width = 10, height = 12, units = "in")
  
  ## Plot 4: Recent events. Smooth lines showing average values for each test statistic as the number of trees increases, faceted by tree number of taxa ##
  # Set dataset for plot
  plot_df <- exp2_long_df
  plot_df <- plot_df[plot_df$recombination_type == "Recent", ]
  # Construct plot
  p <- ggplot(plot_df, aes(x = recombination_value, y = value, color = as.factor(tree_depth))) + 
    geom_smooth() + 
    facet_grid(var_label~num_taxa, scales = "fixed", labeller = label_parsed) +
    scale_x_continuous(name = "Proportion of recombinant DNA", breaks = seq(0,0.5, 0.1), labels = seq(0,0.5, 0.1), minor_breaks = seq(0,0.5, 0.05)) +
    scale_y_continuous(name = "Test statistic value", limits = c(0,1.10), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1)) +
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
  plot_path <- paste0(plot_directory, "exp2_plot4_main.figure_num_taxa_Recent.pdf")
  ggsave(p, filename = plot_path, width = 10, height = 12, units = "in")
}



#### 7. Prepare data from empirical datasets ####
if (plot_empirical == TRUE){
  # Open data file from Experiment 1 as a dataframe
  emp_data_file <- grep("Oaks2011", grep("treelikeness_metrics_collated_results", data_files, value = TRUE), value = TRUE)
  emp_df <- read.csv(file = emp_data_file, stringsAsFactors = FALSE)
  # Convert TIGER values to numeric
  emp_df$mean_TIGER_value <- as.numeric(emp_df$mean_TIGER_value)
  # Convert sCF values to decimal from percentage
  emp_df$sCF_mean <- emp_df$sCF_mean / 100
  
  # Remove rows with ALL codon positions
  emp_df <- emp_df[emp_df$codon_position != "All",]
  
  # Create wide df with only necessary columns
  emp_wide_df <- emp_df[, c("row_id", "uid", "gene_name", "num_taxa", "codon_position", "DNA_type",
                            "tree_proportion", "Cunningham_test", "mean_delta_plot_value", 
                            "LM_proportion_resolved_quartets","NetworkTreelikenessTest",
                            "mean_Q_residual", "sCF_mean", "mean_TIGER_value")]
  
  # Rename the Network Treelikeness Test results
  emp_wide_df$NetworkTreelikenessTest <- emp_wide_df$NetworkTreelikenessTest
  emp_wide_df$NetworkTreelikenessTest[emp_wide_df$NetworkTreelikenessTest == "Zero_splits_where_confidence_intervals_exclude_0"] <- 0
  emp_wide_df$NetworkTreelikenessTest[emp_wide_df$NetworkTreelikenessTest == "Treelike"] <- 1
  emp_wide_df$NetworkTreelikenessTest[emp_wide_df$NetworkTreelikenessTest == "Non-treelike"] <- 0
  
  # Melt exp1_wide_df for better plotting
  emp_long_df <- melt(emp_wide_df, id.vars = c("row_id", "uid", "gene_name", "num_taxa", "codon_position", "DNA_type"))
  emp_long_df$value <- as.numeric(emp_long_df$value)
  
  # Add fancy labels for facets
  emp_long_df$var_label <- factor(emp_long_df$variable, 
                                  levels = c("tree_proportion", "Cunningham_test", "mean_delta_plot_value", 
                                             "LM_proportion_resolved_quartets","NetworkTreelikenessTest",
                                             "mean_Q_residual", "sCF_mean", "mean_TIGER_value"), 
                                  ordered = TRUE, 
                                  labels = c(expression(atop("Tree","proportion")), expression(atop("Cunningham","metric")), 
                                             expression(paste('Mean ', delta["q"])), expression(atop("Proportion","resolved quartets")),
                                             expression(atop("Proportion","treelike alignments")), expression(atop("Mean", "Q-Residual value")), 
                                             expression(atop("Mean", "sCF value")), expression(atop("Mean","TIGER value"))) )
}



#### 8. Plot data from empirical datasets ####
if (plot_empirical == TRUE){
  ## Plot 1: box plots of average value for each test statistic against codon position (faceted by number of taxa, coloured by type of DNA alignment) ## 
  # Set dataset for plot
  plot_df <- emp_long_df

  # Add color blind friendly palette
  colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # grey, orange, sky blue, bluish green, yellow, blue, vermilion, reddish purple
  # Construct plot
  p <- ggplot(plot_df, aes(x = codon_position, y = value, fill = DNA_type)) + 
    geom_boxplot() +
    facet_grid(var_label~num_taxa, scales = "fixed", labeller = label_parsed) +
    scale_x_discrete(name = "Codon position") +
    scale_y_continuous(name = "Test statistic value")+
    guides(fill = guide_legend(title = "Type of DNA\nalignment")) +
    labs(title = "Number of taxa") +
    scale_fill_manual(values = c(colorBlindGrey8[2], colorBlindGrey8[6]), breaks = c("mtDNA", "nDNA"), labels = c("mtDNA", "nDNA")) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14),
          legend.title = element_text(size = 16, hjust = 0.5), legend.text = element_text(size = 14),
          strip.text = element_text(size = 11),
          plot.title = element_text(size = 16, hjust = 0.5))
  # Save plot
  plot_path <- paste0(plot_directory, "emp_plot1_main.figure.pdf")
  ggsave(p, filename = plot_path, width = 10, height = 12, units = "in")
}





