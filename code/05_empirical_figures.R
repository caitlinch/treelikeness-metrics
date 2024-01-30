# caitlinch/treelikeness-metrics/code/05_empirical_figures.R
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
             measure.vars = c("tree_proportion_p_value", "sCF_mean_p_value", "mean_delta_plot_p_value"))



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
                              levels = c("tree_proportion_p_value", "sCF_mean_p_value", "mean_delta_plot_p_value"),
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
quilt_file <- paste0(plot_directory, "mainfig_gene_treelikeness_composite")
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



#### 7. Compare distributions of test statistic values ####
### For tree proportion:
## F test (compare variance)
# F = 1.3279, num df = 116, denom df = 41, p-value = 0.3015
# Accept null hypothesis that variances are equal
tp_f <- var.test(gene_op_df$tree_proportion_alignment_value[which(gene_op_df$dataset == "WEA17")],
                 gene_op_df$tree_proportion_alignment_value[which(gene_op_df$dataset == "WEA17F")], 
                 alternative = "two.sided")
## t-test (compare means)
# t = -3.294, df = 157, p-value = 0.001221
# Reject null hypothesis that the true difference in means is equal to 0
tp_t <- t.test(gene_op_df$tree_proportion_alignment_value[which(gene_op_df$dataset == "WEA17")], 
               gene_op_df$tree_proportion_alignment_value[which(gene_op_df$dataset == "WEA17F")],
               alternative = c("two.sided"),
               paired = FALSE, var.equal = TRUE,
               conf.level = 0.95)

### For mean sCF
## F test (compare variance)
# F = 2.3656, num df = 116, denom df = 41, p-value = 0.002256
# Reject null hypothesis that variances are equal (in t.test function, set var.equal to FALSE)
scf_f <- var.test(gene_op_df$sCF_mean_alignment_value[which(gene_op_df$dataset == "WEA17")],
                  gene_op_df$sCF_mean_alignment_value[which(gene_op_df$dataset == "WEA17F")], 
                  alternative = "two.sided")
## t-test (compare means)
# t = 0.52804, df = 111.72, p-value = 0.5985
# Accept null hypothesis that the true difference in means is equal to 0
scf_t <- t.test(gene_op_df$sCF_mean_alignment_value[which(gene_op_df$dataset == "WEA17")], 
                gene_op_df$sCF_mean_alignment_value[which(gene_op_df$dataset == "WEA17F")],
                alternative = c("two.sided"),
                paired = FALSE, var.equal = FALSE,
                conf.level = 0.95)


### For mean delta plot
## F test (compare variance)
# F = 1.323, num df = 116, denom df = 41, p-value = 0.308
# Accept null hypothesis that variances are equal
dp_f <- var.test(gene_op_df$mean_delta_plot_alignment_value[which(gene_op_df$dataset == "WEA17")],
                 gene_op_df$mean_delta_plot_alignment_value[which(gene_op_df$dataset == "WEA17F")], 
                 alternative = "two.sided")
## t-test (compare means)
# t = 2.9828, df = 157, p-value = 0.003312
# Reject null hypothesis that the true difference in means is equal to 0
dp_t <- t.test(gene_op_df$mean_delta_plot_alignment_value[which(gene_op_df$dataset == "WEA17")], 
               gene_op_df$mean_delta_plot_alignment_value[which(gene_op_df$dataset == "WEA17F")],
               alternative = c("two.sided"),
               paired = FALSE, var.equal = TRUE,
               conf.level = 0.95)



### 8. Plot violin plots for manuscript ####
# Create long format dataframe for p-values
violin_df <- melt(gene_op_df,
                  id.vars = c("dataset", "gene"),
                  measure.vars = c("tree_proportion_alignment_value", "sCF_mean_alignment_value", "mean_delta_plot_alignment_value"))
# Create labels for plotting
violin_df$dataset_label <- factor(violin_df$dataset,
                                  levels = unique(violin_df$dataset),
                                  labels = c("Original", "Orthology-enriched"),
                                  ordered = TRUE)
violin_df$variable_label <- factor(violin_df$variable,
                                   levels = c("tree_proportion_alignment_value", "sCF_mean_alignment_value", "mean_delta_plot_alignment_value"),
                                   labels = c("Tree proportion", "Mean sCF", "Mean delta plot"),
                                   ordered = TRUE)
# Plot the bar chart for the test statistics
violin_plot <- ggplot(data = violin_df, aes(x = variable_label, y = value, fill = dataset_label)) +
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
violin_path <- paste0(plot_directory, "gene_test_statistics_violin")
ggsave(filename = paste0(violin_path, ".png"), plot = violin_plot)
ggsave(filename = paste0(violin_path, ".pdf"), plot = violin_plot)



