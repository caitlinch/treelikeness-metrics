# caitlinch/treelikeness-metrics/code/04_empirical_figures.R
# Caitlin Cherryh 2023

## Script summary:
# This program will plot the treelikeness metrics and p-values for genes in two empirical phylogenetic alignments

#### 1. Set parameters ####
## Directories
# data_directory      <- Directory containing output from gene treelikeness metrics and gene parametric bootstraps
# figure_directory    <- Directory for saving figures
# repo_directory      <- Location of caitlinch/treelikeness-metrics github repository (for access to functions)


# Directories
empirical_data_directory        <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/03_empirical_tree_estimation/"
output_directory                <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/05_empirical_treelikeness_results/"
replicate_alignment_directory   <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/04_bs_replicate_alignments/"
repo_directory                  <- "/Users/caitlincherryh/Documents/Repositories/treelikeness-metrics/"



#### 2. Prepare analyses ####
# Open packages
library(reshape2)
library(ggplot2)
library(patchwork)
library(viridisLite)

# Source functions from caitlinch/treelikeness_metrics
source(paste0(repo_directory, "code/func_empirical.R"))



#### 6. Plot histograms of empirical bootstrap test statistic values ####
if (control_parameters$plots == TRUE){
  ## Find file with parametric bootstrap results
  all_output <- list.files(output_directory, recursive = TRUE)
  all_csv <- grep("csv", all_output, value = T)
  bs_csv_file <- paste0(output_directory, grep("empirical", grep("collated", grep("treelikeness_metrics", all_csv, value = T), value = T), value = T))
  # Open csv file
  bs_df <- read.csv(bs_csv_file, stringsAsFactors = F)
  
  #### Plot 1: histogram of test statistic values ####
  ## Reformat dataframe for plotting
  # Add a column specifying whether the row is for an alignment or a bootstrap replicate
  bs_df$rep_type <- factor(bs_df$unique_id,
                           levels = bs_df$unique_id,
                           labels = rep(c(rep("Bootstrap replicate", 100), "Alignment"), 2),
                           ordered = TRUE)
  # Add a column specifying which alignment the row is
  bs_df$alignment_id <- factor(bs_df$unique_id,
                               levels = bs_df$unique_id,
                               labels = rep(c("WEA17", "WEA17F"), each = 101),
                               ordered = TRUE)
  # Divide scf_mean by 100 to make it a proportion
  bs_df$sCF_mean <- bs_df$sCF_mean/100
  
  # # Copy dataframe and remove alignment rows
  # bs_only_df <- bs_df[which(bs_df$rep_type == "Bootstrap replicate"), ]
  # 
  # # Reformat dataframe into long format
  # long_df <- melt(bs_only_df,
  #                 id.vars = c("alignment_id", "rep_type", "unique_id"),
  #                 measure.vars = c("LM_proportion_resolved_quartets", "sCF_mean", "mean_delta_plot_value", "Cunningham_test", "tree_proportion"))
  
  # Reformat dataframe into long format
  long_df <- melt(bs_df,
                  id.vars = c("alignment_id", "rep_type", "unique_id"),
                  measure.vars = c("LM_proportion_resolved_quartets", "sCF_mean", "mean_delta_plot_value", "Cunningham_test", "tree_proportion"))
  # Replace all Alignment values with NA to avoid plotting them in histogram - want to plot them as vertical line instead
  long_df$value[which(long_df$rep_type == "Alignment")] <- NA
  
  ## Add variables for facetting labels nicely
  long_df$var_label = factor(long_df$variable,
                             levels = c("tree_proportion", "Cunningham_test", "mean_delta_plot_value", "LM_proportion_resolved_quartets", "sCF_mean"),
                             ordered = TRUE,
                             labels = c(expression(atop("Tree","proportion")), 
                                        expression(atop("Cunningham","metric")),
                                        expression(paste('Mean ', delta["q"])), 
                                        expression(atop(textstyle("Proportion"),atop(textstyle("resolved"),atop(scriptscriptstyle(""),textstyle("quartets"))))),
                                        expression(atop("Mean", "sCFL value"))) )
  long_df$alignment_label <- factor(long_df$alignment_id,
                                    levels = c("WEA17", "WEA17F"),
                                    ordered = TRUE,
                                    labels = c("Whelan 2017\nOriginal dataset", "Filtered by\nMcCarthy 2023") )
  
  ## Add variable detailing each test and the empirical value for each test
  # Add new row for combination of alignment and test
  xinterval_df <- data.frame(alignment_id = c(rep(c("WEA17",  "WEA17F"), each = 5)),
                             rep_type = rep("Alignment", 10),
                             unique_id = c(rep(c("WEA17",  "WEA17F"), each = 5)),
                             variable = rep(c("tree_proportion", "Cunningham_test", "mean_delta_plot_value", "LM_proportion_resolved_quartets", "sCF_mean"), 2),
                             value = c(bs_df[bs_df$alignment_id == "WEA17" & bs_df$rep_type == "Alignment",]$tree_proportion,
                                       bs_df[bs_df$alignment_id == "WEA17" & bs_df$rep_type == "Alignment",]$Cunningham_test,
                                       bs_df[bs_df$alignment_id == "WEA17" & bs_df$rep_type == "Alignment",]$mean_delta_plot_value,
                                       bs_df[bs_df$alignment_id == "WEA17" & bs_df$rep_type == "Alignment",]$LM_proportion_resolved_quartets,
                                       bs_df[bs_df$alignment_id == "WEA17" & bs_df$rep_type == "Alignment",]$sCF_mean,
                                       bs_df[bs_df$alignment_id == "WEA17F" & bs_df$rep_type == "Alignment",]$tree_proportion,
                                       bs_df[bs_df$alignment_id == "WEA17F" & bs_df$rep_type == "Alignment",]$Cunningham_test,
                                       bs_df[bs_df$alignment_id == "WEA17F" & bs_df$rep_type == "Alignment",]$mean_delta_plot_value,
                                       bs_df[bs_df$alignment_id == "WEA17F" & bs_df$rep_type == "Alignment",]$LM_proportion_resolved_quartets,
                                       bs_df[bs_df$alignment_id == "WEA17F" & bs_df$rep_type == "Alignment",]$sCF_mean))
  xinterval_df$var_label = factor(xinterval_df$variable,
                                  levels = c("tree_proportion", "Cunningham_test", "mean_delta_plot_value", "LM_proportion_resolved_quartets", "sCF_mean"),
                                  ordered = TRUE,
                                  labels = c(expression(atop("Tree","proportion")), 
                                             expression(atop("Cunningham","metric")),
                                             expression(paste('Mean ', delta["q"])), 
                                             expression(atop(textstyle("Proportion"),atop(textstyle("resolved"),atop(scriptscriptstyle(""),textstyle("quartets"))))),
                                             expression(atop("Mean", "sCFL value"))) )
  xinterval_df$alignment_label <- factor(xinterval_df$alignment_id,
                                         levels = c("WEA17", "WEA17F"),
                                         ordered = TRUE,
                                         labels = c("Whelan 2017\nOriginal dataset", "Filtered by\nMcCarthy 2023") )
  
  ## Set colour palette 
  #   To use scale_fill_manual:     scale_fill_manual(labels = c("Bootstrap replicate", "Alignment"), values = c(col_pal[2], col_pal[1]))
  col_pal <- cividis(2, begin = 0, end = 1, direction = 1) # Note: col_pal = c("#00204DFF", "#FFEA46FF") = c("dark blue", "bright yellow")
  
  ## Add the expression for the plot title
  hist_title_expression <- expression(atop("Test statistic values for Whelan et al 2017","alignments and parametric bootstrap replicates"))
  ## Plot a nice histogram of the output values
  h <- ggplot(long_df, aes(x = value, fill = rep_type, colour=rep_type)) +
    geom_histogram(bins = 20) +
    geom_vline(data = xinterval_df, aes(xintercept = value), colour = "#00204DFF", show.legend = F, linetype = "dashed") + 
    facet_grid(alignment_label~var_label, labeller = labeller(var_label = label_parsed)) +
    scale_fill_manual(labels = c("Bootstrap replicate", "Alignment"),
                      values = c(col_pal[2], col_pal[1]),
                      name = "Legend") +
    scale_colour_manual(labels = c("Bootstrap replicate", "Alignment"),
                        values = c(col_pal[2], col_pal[1]),
                        name = "Legend") +  
    ggtitle(parse(text = hist_title_expression)) +
    scale_x_continuous(name = "Test statistic value") +
    scale_y_continuous(name = "Count") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title.x = element_text(size = 18, margin = margin(t = 15, r = 0, b = 10, l = 0)),
          axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 15, b = 0, l = 10)), 
          axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = 12),
          strip.text.x = element_text(size = 12),
          strip.text.y = element_text(size = 12),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 12)) 
  ## Save the histogram
  h_path <- paste0(output_directory, "empirical_test_statistic_value_histogram")
  ggsave(filename = paste0(h_path, ".png"), plot = h)
  ggsave(filename = paste0(h_path, ".pdf"), plot = h)
  
  #### Plot 2: proportion of treelike alignments ####
  # Make a new dataframe with the proportion of treelike alignments
  prop_tl_alignments <- c(length(which(bs_df$rep_type == "Bootstrap replicate" & bs_df$alignment_id == "WEA17" & bs_df$NetworkTreelikenessTest == "Treelike")),
                          length(which(bs_df$rep_type == "Bootstrap replicate" & bs_df$alignment_id == "WEA17" & bs_df$NetworkTreelikenessTest == "Non-treelike")),
                          length(which(bs_df$rep_type == "Bootstrap replicate" & bs_df$alignment_id == "WEA17F" & bs_df$NetworkTreelikenessTest == "Treelike")),
                          length(which(bs_df$rep_type == "Bootstrap replicate" & bs_df$alignment_id == "WEA17F" & bs_df$NetworkTreelikenessTest == "Non-treelike")) )
  tl_df <- data.frame(alignment_id = c(rep("WEA17", prop_tl_alignments[1]), rep("WEA17", prop_tl_alignments[2]),rep("WEA17F", prop_tl_alignments[3]), rep("WEA17F", prop_tl_alignments[4]), "WEA17", "WEA17F"),
                      rep_type = c(rep("Bootstrap replicate", prop_tl_alignments[1]), rep("Bootstrap replicate", prop_tl_alignments[2]),rep("Bootstrap replicate", prop_tl_alignments[3]), rep("Bootstrap replicate", prop_tl_alignments[4]), "Alignment", "Alignment"),
                      test_value = c(rep("Treelike", prop_tl_alignments[1]), rep("Non-treelike", prop_tl_alignments[2]),rep("Treelike", prop_tl_alignments[3]), rep("Non-treelike", prop_tl_alignments[4]), "Non-treelike", "Non-treelike") )
  tl_df$test_value <- factor(tl_df$test_value,
                             levels = c("Treelike", "Non-treelike", "Treelike", "Non-treelike", "Non-treelike", "Non-treelike"),
                             ordered = T,
                             labels = c("Treelike", "Non-treelike", "Treelike", "Non-treelike", "Non-treelike", "Non-treelike") )
  tl_df$alignment_label <- factor(tl_df$alignment_id,
                                  levels = c("WEA17", "WEA17F"),
                                  ordered = TRUE,
                                  labels = c("Whelan 2017\nOriginal dataset", "Filtered by\nMcCarthy 2023") )
  # Remove alignment rows
  tl_df <- tl_df[tl_df$rep_type == "Bootstrap replicate", ]
  
  ## Add the expression for the plot title
  bar_title_expression <- expression(atop("Network Treelikeness Test values for Whelan et al 2017","alignments and parametric bootstrap replicates"))
  ## Plot a nice histogram of the output values
  b <- ggplot(tl_df, aes(x = test_value, fill = rep_type)) +
    geom_bar() +
    facet_grid(~alignment_label) +
    scale_fill_viridis_d(option = "E", direction = -1) +
    guides(fill=guide_legend(title="Legend")) +
    ggtitle(parse(text = bar_title_expression)) +
    scale_x_discrete(name = "Value") +
    scale_y_continuous(name = "Count") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 18),
          axis.title.x = element_text(size = 18, margin = margin(t = 15, r = 0, b = 10, l = 0)),
          axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 15, b = 0, l = 10)), 
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 12),
          strip.text.x = element_text(size = 12),
          strip.text.y = element_text(size = 12),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 12))
  ## Save the bar plot
  b_path <- paste0(output_directory, "empirical_test_statistic_value_bar")
  ggsave(filename = paste0(b_path, ".png"), plot = b)
  ggsave(filename = paste0(b_path, ".pdf"), plot = b)
  
  #### Plot 3: combine plots 1 and 2 with patchwork ####
  ## Collate the two plots using patchwork
  quilt = h / b + plot_annotation(tag_levels = 'a', tag_suffix = '.') & 
    theme(plot.tag = element_text(size = 20))
  ## Save the collated plot
  quilt_path <- paste0(output_directory, "empirical_collated_plot")
  ggsave(filename = paste0(quilt_path, ".png"), plot = quilt, width = 9, height = 12, units = "in")
  ggsave(filename = paste0(quilt_path, ".pdf"), plot = quilt, width = 9, height = 12, units = "in")
}
