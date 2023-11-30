# caitlinch/treelikeness-metrics/code/04_empirical_figures.R
# Caitlin Cherryh 2023

## Script summary:
# This program will plot the treelikeness metrics and p-values for genes in two empirical phylogenetic alignments

#### 1. Set parameters ####
## Directories
# figure_directory    <- Directory for saving figures
# repo_directory      <- Location of caitlinch/treelikeness-metrics github repository (for access to functions)

# Directories
repo_directory <- "/Users/caitlincherryh/Documents/Repositories/treelikeness-metrics/"
plot_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/06_plots/"



#### 2. Prepare analyses ####
# Open packages
library(reshape2)
library(ggplot2)
library(patchwork)

# Open empirical gene p-value file
all_op_files <- paste0(paste0(repo_directory, "output/"), list.files(paste0(repo_directory, "output/")))
gene_op_file <- grep("pvalues", all_op_files, value = T)
# Open gene output df
gene_op_df <- read.csv(gene_op_file, stringsAsFactors = FALSE)
# Scale sCF_mean_alignment_value and sCF_median_alignment_value to match other test statistics range of 0-1
gene_op_df$sCF_mean_alignment_value <- gene_op_df$sCF_mean_alignment_value/100
gene_op_df$sCF_median_alignment_value <- gene_op_df$sCF_median_alignment_value/100
# Create long format dataframe for test statistics
ts_df <- melt(gene_op_df,
              id.vars = c("dataset", "gene"),
              measure.vars = c("tree_proportion_alignment_value", "sCF_mean_alignment_value", "mean_delta_plot_alignment_value"))
# Create long format dataframe for p-values
p_df <- melt(gene_op_df,
             id.vars = c("dataset", "gene"),
             measure.vars = c("tree_proportion_p_value_ecdf", "sCF_mean_p_value_ecdf", "mean_delta_plot_p_value_ecdf"))



#### 3. Plot test statistic histograms ####
# Create new column for whether p-value <= 0.05
ts_df$p_value_significance <- p_df$value
ts_df$p_value_significance[p_df$value <= 0.05] <- "Significant"
ts_df$p_value_significance[p_df$value > 0.05] <- "Not significant"
# Create labels for plotting
ts_df$dataset_label <- factor(ts_df$dataset,
                              levels = unique(ts_df$dataset),
                              labels = c("Original dataset", "Orthology-enriched\ndataset"),
                              ordered = TRUE)
ts_df$variable_label <- factor(ts_df$variable,
                               levels = c("tree_proportion_alignment_value", "sCF_mean_alignment_value", "mean_delta_plot_alignment_value"),
                               labels = c("Tree proportion", "Mean sCF", "Mean delta plot"),
                               ordered = TRUE)
# Plot histogram for each test statistic, faceted by test statistic and dataset
hist <- ggplot(data = ts_df, aes(x = value)) + 
  geom_histogram(binwidth = 0.05) +
  facet_grid(dataset_label ~ variable_label, scale = "fixed") +
  scale_y_continuous(name = "Frequency") +
  scale_x_continuous(name = "Test statistic value", limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), minor_breaks = seq(0, 1, 0.1)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 10, l = 0, unit = "pt")),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 0.5, vjust = 0.5), 
        axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 10, b = 0, l = 10, unit = "pt")), 
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 13, margin = margin(t = 5, r = 0, b = 5, l = 0, unit = "pt")),
        strip.text.y = element_text(size = 13, margin = margin(t = 0, r = 5, b = 0, l = 5, unit = "pt")))



#### 4. Plot p-value histograms ####
# Create new column for whether p-value <= 0.05
p_df$p_value_significance <- p_df$value
p_df$p_value_significance[p_df$value <= 0.05] <- "Significant"
p_df$p_value_significance[p_df$value > 0.05] <- "Not significant"
# Create labels for plotting
p_df$dataset_label <- factor(p_df$dataset,
                             levels = unique(p_df$dataset),
                             labels = c("Original dataset", "Orthology-enriched dataset"),
                             ordered = TRUE)
p_df$variable_label <- factor(p_df$variable,
                              levels = c("tree_proportion_p_value_ecdf", "sCF_mean_p_value_ecdf", "mean_delta_plot_p_value_ecdf"),
                              labels = c("Tree proportion", "Mean sCF", "Mean delta plot"),
                              ordered = TRUE)
# Plot histogram for each p-value, faceted by test statistic and dataset
p_hist <- ggplot(data = p_df, aes(x = value)) + 
  geom_histogram(bins = 20) +
  facet_grid(variable_label ~ dataset_label, scale = "fixed") +
  scale_x_continuous(name = "p-value", breaks = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), minor_breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(name = "Frequency", breaks = seq(0, 120, 30), labels = seq(0, 120, 30), minor_breaks = seq(0, 120, 15), limits = c(0, 120)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 10, l = 0, unit = "pt")),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 0.5, vjust = 0.5), 
        axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 10, b = 0, l = 10, unit = "pt")), 
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 16, margin = margin(t = 5, r = 0, b = 5, l = 0, unit = "pt")),
        strip.text.y = element_text(size = 16, margin = margin(t = 0, r = 5, b = 0, l = 5, unit = "pt")))



#### 5. Plot p-value stacked bar charts ####
# Plot stacked bar charts
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

# Create tables
table(og_p_df$p_value_significance, og_p_df$variable_label)
table(filtered_p_df$p_value_significance, filtered_p_df$variable_label)



#### 6. Assemble plots using patchwork ####
# Assemble main figure
quilt <- hist/(og_p_plot + filtered_p_plot) + 
  plot_annotation(tag_levels = 'a', tag_suffix = ")") &  theme(plot.tag = element_text(size = 25))
quilt_file <- paste0(plot_directory, "gene_treelikeness_composite")
ggsave(filename = paste0(quilt_file, ".pdf"), plot = quilt, width = 10, height = 12)
ggsave(filename = paste0(quilt_file, ".png"), plot = quilt, width = 10, height = 12)
# Save test statistic histogram
plot_file <- paste0(plot_directory, "gene_treelikeness_test_statistic_histogram")
ggsave(filename = paste0(plot_file, ".pdf"), plot = hist)
ggsave(filename = paste0(plot_file, ".png"), plot = hist)
# Save p-value histogram
plot_file <- paste0(plot_directory, "gene_treelikeness_pvalue_histogram")
ggsave(filename = paste0(plot_file, ".pdf"), plot = p_hist)
ggsave(filename = paste0(plot_file, ".png"), plot = p_hist)
# Save stacked box plots of p-values
quilt2 <- og_p_plot + filtered_p_plot + plot_layout(ncol = 2, nrow = 1) + 
  plot_annotation(tag_levels = 'a', tag_suffix = ")") &  theme(plot.tag = element_text(size = 25))
quilt2_file <- paste0(plot_directory, "gene_treelikeness_pvalue_stacked_bar")
ggsave(filename = paste0(quilt2_file, ".pdf"), plot = quilt2)
ggsave(filename = paste0(quilt2_file, ".png"), plot = quilt2)


