# caitlinch/treelikeness-metrics/code/06_presentation_figures.R
# Caitlin Cherryh 2023


## Script summary:
# This program will plot the treelikeness metrics and p-values for genes in two empirical phylogenetic alignments

#### 1. Set parameters ####
## Directories
# repo_directory      <- Location of caitlinch/treelikeness-metrics github repository (for access to functions)
# plot_directory      <- Directory for saving figures

# Directories
repo_directory <- "/Users/caitlincherryh/Documents/Repositories/treelikeness-metrics/"
plot_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/06_plots/"



#### 2. Prepare analyses ####
# Open packages
library(reshape2)
library(ggplot2)
library(patchwork)

# List all output files
output_files <- paste0(paste0(repo_directory, "output/"), list.files(paste0(repo_directory, "output/")))



### 3. Plot tree proportion for random tree simulations ####
# Open random tree simulations
rt_file_path <- grep("collated_results", grep("randomTree", output_files, value = T), value = T)
rt_df <- read.csv(rt_file_path)
# Remove columns you don't want for plotting
rt_wide_df <- rt_df[, c("row_id", "uid", "num_taxa", "num_trees", "tree_depth", 
                        "tree_proportion")]
# Melt exp1_wide_df for better plotting
rt_long_df <- melt(rt_wide_df, id.vars = c("row_id", "uid", "num_taxa", "num_trees", "tree_depth"), measure.vars = c("tree_proportion"))
# Add fancy labels for facets
rt_long_df$var_label <- factor(rt_long_df$variable, 
                               levels = c("tree_proportion"), 
                               ordered = TRUE, 
                               labels = c(expression(atop("Tree","proportion"))) )
# Set log10 minor breaks for x and y axis
x_axis_minor_breaks <-  unique(c(seq(1, 10, 1), seq(10, 100, 10), seq(100, 1000, 100), seq(1000, 10000, 1000)))
# Create plot for presentation
p1 <- ggplot(rt_long_df, aes(x = num_trees, y = value, color = as.factor(num_taxa))) + 
  geom_smooth(method = "loess", alpha = 0.2, linewidth = 0, span = 0.75,
              aes(ymin = ifelse(after_stat(ymin) < 0, 0, after_stat(ymin)), ymax = ifelse(after_stat(ymax) > 1, 1, after_stat(ymax)))) +
  stat_smooth(method = "loess", geom = "line", linewidth = 1.2, alpha = 0.7, span = 0.75) +
  facet_grid(var_label~tree_depth, scales = "fixed", labeller = label_parsed) +
  scale_x_log10(minor_breaks = x_axis_minor_breaks) +
  scale_y_continuous(name = "Test statistic value", limits = c(0,1.10), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1), oob=scales::rescale_none) +
  scale_color_viridis_d(direction = -1) +
  guides(color = guide_legend(title = "Number\nof taxa")) +
  labs(subtitle = "Tree depth (substitutions per site)", 
       x = expression("Number of trees ("*log[10]*" scale)")) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 22, margin = margin(t = 20, r = 0, b = 10, l = 0), color = "grey30"), 
        axis.text.x = element_text(size = 18, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 22, margin = margin(t = 0, r = 20, b = 0, l = 10), color = "grey30"), 
        axis.text.y = element_text(size = 18),
        legend.title = element_text(size = 22), 
        legend.text = element_text(size = 18),
        strip.text.y = element_blank(),
        strip.text.x = element_text(size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 22, margin = margin(t = 10, r = 0, b = 20, l = 0)) )
# Save the plot
p1_path <- paste0(plot_directory, "es_tree_proportion_randomTrees")
ggsave(filename = paste0(p1_path, ".png"), plot = p1)
ggsave(filename = paste0(p1_path, ".pdf"), plot = p1)



### 4. Plot empirical values (tree proportion) ####
# Open gene empirical results
gene_file_path <- grep("pvalues", grep("gene", output_files, value = T), value = T)
gene_df <- read.csv(gene_file_path)
# Scale sCF values
gene_df$sCF_mean_alignment_value <- gene_df$sCF_mean_alignment_value/100
# Create long format dataframe for p-values
ts_df <- melt(gene_df,
             id.vars = c("dataset", "gene"),
             measure.vars = c("tree_proportion_alignment_value", "sCF_mean_alignment_value", "mean_delta_plot_alignment_value"))
# Create labels for plotting
ts_df$dataset_label <- factor(ts_df$dataset,
                             levels = unique(ts_df$dataset),
                             labels = c("Original", "Orthology-enriched"),
                             ordered = TRUE)
ts_df$variable_label <- factor(ts_df$variable,
                              levels = c("tree_proportion_alignment_value", "sCF_mean_alignment_value", "mean_delta_plot_alignment_value"),
                              labels = c("Tree proportion", "Mean sCF", "Mean delta plot"),
                              ordered = TRUE)
# Plot the bar chart for the test statistics
p2 <- ggplot(data = ts_df, aes(x = variable_label, y = value, fill = dataset_label)) +
  geom_boxplot() +
  scale_fill_manual(name = "Dataset", values = c("#d01c8b", "#b8e186")) +
  scale_y_continuous(name = "Test statistic value", limits = c(0, 1), breaks = seq(0,1,0.2), labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.05)) +
  scale_x_discrete(name = NULL) +
  guides(fill = guide_legend(override.aes = list(size = 12))) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20, angle = 45, vjust = 1.0, hjust = 1.0, margin = margin(t = 5, r = 0, b = 0, l = 0, unit = "pt")), 
        axis.title.y = element_text(size = 26, margin = margin(t = 0, r = 20, b = 0, l = 10, unit = "pt"), color = "grey30"), 
        axis.text.y = element_text(size = 20, margin = margin(t = 0, r = 5, b = 0, l = 0, unit = "pt")),
        legend.title = element_text(size = 26),
        legend.text = element_text(size = 24),
        axis.line = element_line(linewidth = 1), 
        axis.ticks = element_line(linewidth = 1),
        panel.grid.major = element_line(linewidth = 1),
        panel.grid.minor = element_line(linewidth = 0.8),
        panel.border = element_rect(linewidth = 2, fill = NA))
# Save the plot
p2_path <- paste0(plot_directory, "es_gene_test_statistics")
ggsave(filename = paste0(p2_path, ".png"), plot = p2)
ggsave(filename = paste0(p2_path, ".pdf"), plot = p2)



### 5. Plot empirical values (p-values) ####
# Open gene empirical results
gene_file_path <- grep("pvalues", grep("gene", output_files, value = T), value = T)
gene_df <- read.csv(gene_file_path)
# Create long format dataframe for p-values
p_df <- melt(gene_df,
             id.vars = c("dataset", "gene"),
             measure.vars = c("tree_proportion_p_value_ecdf", "sCF_mean_p_value_ecdf", "mean_delta_plot_p_value_ecdf"))
# Create new column for whether p-value <= 0.05
p_df$p_value_significance <- p_df$value
p_df$p_value_significance[p_df$value <= 0.05] <- "Significant"
p_df$p_value_significance[p_df$value > 0.05] <- "Not significant"
# Create labels for plotting
p_df$dataset_label <- factor(p_df$dataset,
                             levels = unique(p_df$dataset),
                             labels = c("Original dataset", "Orthology-enriched\ndataset"),
                             ordered = TRUE)
p_df$variable_label <- factor(p_df$variable,
                              levels = c("tree_proportion_p_value_ecdf", "sCF_mean_p_value_ecdf", "mean_delta_plot_p_value_ecdf"),
                              labels = c("Tree proportion", "Mean sCF", "Mean delta plot"),
                              ordered = TRUE)
# Separate dataframes for separate datasets
og_p_df <- p_df[which(p_df$dataset == "WEA17"), ]
filtered_p_df <- p_df[which(p_df$dataset == "WEA17F"), ]
# Plot original dataset
og_p_plot <- ggplot(data = og_p_df, aes(x = variable_label, fill = p_value_significance)) +
  geom_bar() +
  scale_fill_manual(name = "p-value result", values = c("#67a9cf", "#ef8a62")) +
  labs(title = "Original dataset") +
  scale_x_discrete(name = "Test statistic p-values") +
  scale_y_continuous(name = "Count") +
  guides(fill = "none") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16, margin = margin(t = 5, r = 0, b = 5, l = 0, unit = "pt")),
        axis.text.x = element_text(size = 14, margin = margin(t = 5, r = 0, b = 5, l = 0, unit = "pt"), hjust = 1, vjust = 1, angle = 45), 
        axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 10, b = 0, l = 10, unit = "pt")), 
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 20, hjust = 0.5, margin = margin(t = 10, r = 0, b = 10, l = 0, unit = "pt")) )
# Plot filtered dataset
filtered_p_plot <- ggplot(data = filtered_p_df, aes(x = variable_label, fill = p_value_significance)) +
  geom_bar() +
  scale_fill_manual(name = "p-value result", values = c("#67a9cf", "#ef8a62")) +
  labs(title = "Orthology-enriched dataset") +
  scale_x_discrete(name = "Test statistic p-values") +
  scale_y_continuous(name = "Count") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16, margin = margin(t = 5, r = 0, b = 5, l = 0, unit = "pt")),
        axis.text.x = element_text(size = 14, margin = margin(t = 5, r = 0, b = 5, l = 0, unit = "pt"), hjust = 1, vjust = 1, angle = 45), 
        axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 10, b = 0, l = 10, unit = "pt")), 
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 20, hjust = 0.5, margin = margin(t = 10, r = 0, b = 10, l = 0, unit = "pt")) )
# Collate
quilt <- (og_p_plot + filtered_p_plot)
# Save the plot
p3_path <- paste0(plot_directory, "es_gene_stacked_p_values")
ggsave(filename = paste0(p3_path, ".png"), plot = quilt)
ggsave(filename = paste0(p3_path, ".pdf"), plot = quilt)





