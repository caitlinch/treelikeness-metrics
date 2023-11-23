# caitlinch/treelikeness-metrics/code/04_empirical_data_test.R
# Caitlin Cherryh 2023

## Script summary:
# This program will apply various tests for treelikeness to empirical alignments
# This program requires IQ-Tree2 (2.2-beta or above), fast TIGER, phylogemetric, and SplitsTree (4.17.2 or above).

#### 1. Set parameters ####
## Directories
# output_directory                <- Directory where alignments will be saved/treelikeness metrics will be run
# replicate_alignment_directory   <- Location of bootstrap replicate alignments, generated using Alisim in IQ-Tree
# repo_directory                  <- Location of caitlinch/treelikeness-metrics github repository (for access to functions)

## Executable paths
# iqtree2_path            <- Path to IQ-Tree (requires version 2.2-beta or above)
# splitstree_path         <- Path to SplitsTree 4 version 4.17.2 or above

## Run parameters
# num_cores               <- Number of parallel threads to use at once

run_location = "local"
if (run_location == "local"){
  # Directories
  empirical_data_directory        <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/03_empirical_tree_estimation/"
  output_directory                <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/05_empirical_treelikeness_results/"
  replicate_alignment_directory   <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/04_bs_replicate_alignments/"
  repo_directory                  <- "/Users/caitlincherryh/Documents/Repositories/treelikeness-metrics/"
  
  # Executable paths
  iqtree2_path      <- "iqtree2"
  splitstree_path   <- "/Applications/SplitsTree/SplitsTree.app/Contents/MacOS/JavaApplicationStub"
  
  # Run parameters
  num_cores <- 1
} else if (run_location == "soma"){
  # Directories
  output_directory                <- "/data/caitlin/treelikeness_metrics/empirical_treelikeness_output/"
  replicate_alignment_directory   <- "/data/caitlin/treelikeness_metrics/empirical_alignments/"
  repo_directory                  <- "/data/caitlin/treelikeness_metrics/"
  
  # Executable paths
  iqtree2_path      <- "/data/caitlin/executables/iqtree-2.2.2-Linux/bin/iqtree2"
  splitstree_path   <- "/home/caitlin/splitstree4/SplitsTree"
  
  # Run parameters
  num_cores <- 30
}

control_parameters <- list(edit.replicate.alignments = FALSE,
                           apply.treelikeness.tests = FALSE,
                           calculate.empirical.p_values = FALSE,
                           plots = FALSE,
                           per.gene.analysis = TRUE)



#### 2. Prepare analyses ####
# Open packages
library(parallel)
if (control_parameters$plots == TRUE){
  # Open libraries for plotting, if plotting will be performed
  library(reshape2)
  library(ggplot2)
  library(patchwork)
  library(viridisLite)
}

# Source functions from caitlinch/treelikeness_metrics
source(paste0(repo_directory, "code/func_empirical.R"))
source(paste0(repo_directory, "code/func_parametric_bootstrap.R"))
source(paste0(repo_directory, "code/func_metrics.R"))
source(paste0(repo_directory, "code/func_data_analysis.R"))

if (control_parameters$apply.treelikeness.tests == TRUE){
  # Calculate the number of parallel processes to run using mc.apply
  mclapply_num_cores <- num_cores/10
  if (mclapply_num_cores < 1){
    # Minimum of one process running at once
    mclapply_num_cores <- 1
  } else {
    # Take floor of the decimal to run a conservative number of processes 
    #     (e.g. if 31 cores and each IQ-Tree run uses 10, you can run 3 alignments at once)
    mclapply_num_cores <- floor(mclapply_num_cores)
  }
}



#### 3. Add gaps to the replicate alignments ####
if (control_parameters$edit.replicate.alignments == TRUE){
  ## Extract list of alignments
  all_files <- list.files(replicate_alignment_directory, recursive = T)
  
  ## For WEA17
  # Extract the list of alignments using the id (WEA17F)
  wea17_files <- paste0(replicate_alignment_directory, grep("WEA17_", all_files, value = T))
  wea17_al <- wea17_files[grep("bs_rep", basename(wea17_files), ignore.case = T, invert = T)]
  wea17_replicates <- wea17_files[grep("bs_rep", basename(wea17_files), ignore.case = T)]
  # Apply the gaps from the empirical alignment to the replicate alignments
  lapply(wea17_replicates, copy.alignment.gaps, template_alignment_path = wea17_al)
  
  ## For WEA17_filtered
  # Extract the list of alignments using the id (WEA17F)
  wea17f_files <- paste0(replicate_alignment_directory, grep("WEA17F_", all_files, value = T))
  wea17f_al <- wea17f_files[grep("bs_rep", basename(wea17f_files), ignore.case = T, invert = T)]
  wea17f_replicates <- wea17f_files[grep("bs_rep", basename(wea17f_files), ignore.case = T)]
  # Apply the gaps from the empirical alignment to the replicate alignments
  lapply(wea17f_replicates, copy.alignment.gaps, template_alignment_path = wea17f_al)
}



#### 4. Apply tests for treelikeness to each alignment ####
if (control_parameters$apply.treelikeness.tests == TRUE){
  ## Extract list of alignments
  all_files <- list.files(replicate_alignment_directory, recursive = T)
  
  ## For WEA17
  # Extract the list of alignments using the id (WEA17F)
  wea17_files <- paste0(replicate_alignment_directory, grep("WEA17_", all_files, value = T))
  wea17_al <- wea17_files[grep("bs_rep", basename(wea17_files), ignore.case = T, invert = T)]
  wea17_replicates <- wea17_files[grep("bs_rep", basename(wea17_files), ignore.case = T)]
  # Copy original alignment and prepare for treelikeness metric runs
  wea17_dir <- paste0(output_directory, "WEA17/")
  if (dir.exists(wea17_dir) == FALSE){ dir.create(wea17_dir) }
  wea17_copy <- paste0(wea17_dir, basename(wea17_al))
  file.copy(from = wea17_al, to = wea17_copy)
  # Calculate the tree proportion for the original alignments
  treelikeness.metrics.empirical(wea17_copy,
                                 splitstree_path = splitstree_path, 
                                 iqtree2_path = iqtree2_path, 
                                 iqtree_model = "'Q.insect+R8{0.1639,0.0315,0.1812,0.1846,0.1454,0.4618,0.1171,0.7116,0.1742,1.1628,0.1403,2.0884,0.0676,3.6840,0.0103,6.4538}'", 
                                 num_iqtree2_threads = num_cores, 
                                 sequence_format = "AA", 
                                 redo = FALSE)
  # Apply the wrapper function to calculate treelikeness of bootstrap replicates
  wea17_df <- bootstrap.wrapper(wea17_replicates, 
                                output_directory = output_directory, 
                                splitstree_path = splitstree_path, 
                                iqtree2_path = iqtree2_path, 
                                iqtree_model = "'Q.insect+R8{0.1639,0.0315,0.1812,0.1846,0.1454,0.4618,0.1171,0.7116,0.1742,1.1628,0.1403,2.0884,0.0676,3.6840,0.0103,6.4538}'", 
                                num_iqtree2_threads = "10", 
                                sequence_format = "AA", 
                                redo = FALSE, 
                                number_parallel_cores = mclapply_num_cores)
  
  ## For WEA17_filtered
  # Extract the list of alignments using the id (WEA17F)
  wea17f_files <- paste0(replicate_alignment_directory, grep("WEA17F_", all_files, value = T))
  wea17f_al <- wea17f_files[grep("bs_rep", basename(wea17f_files), ignore.case = T, invert = T)]
  wea17f_replicates <- wea17f_files[grep("bs_rep", basename(wea17f_files), ignore.case = T)]
  # Copy original alignment and prepare for treelikeness metric runs
  wea17f_dir <- paste0(output_directory, "WEA17F/")
  if (dir.exists(wea17f_dir) == FALSE){ dir.create(wea17f_dir) }
  wea17f_copy <- paste0(wea17f_dir, basename(wea17f_al))
  file.copy(from = wea17f_al, to = wea17f_copy)
  # Calculate the tree proportion for the original alignments
  treelikeness.metrics.empirical(wea17f_copy,
                                 splitstree_path = splitstree_path, 
                                 iqtree2_path = iqtree2_path, 
                                 iqtree_model = "'Q.insect+I{0.0659}+R6{0.2004,0.1133,0.1734,0.3827,0.1439,0.6654,0.2108,1.1677,0.1549,2.2361,0.0508,4.3870}'", 
                                 num_iqtree2_threads = num_cores, 
                                 sequence_format = "AA", 
                                 redo = FALSE)
  
  # Apply the wrapper function
  wea17f_df <- bootstrap.wrapper(wea17f_replicates, 
                                 output_directory = output_directory, 
                                 splitstree_path = splitstree_path, 
                                 iqtree2_path = iqtree2_path, 
                                 iqtree_model = "'Q.insect+I{0.0659}+R6{0.2004,0.1133,0.1734,0.3827,0.1439,0.6654,0.2108,1.1677,0.1549,2.2361,0.0508,4.3870}'", 
                                 num_iqtree2_threads = "10", 
                                 sequence_format = "AA", 
                                 redo = FALSE, 
                                 number_parallel_cores = mclapply_num_cores)
}



#### 5. Calculate p-values for each test statistic ####
if (control_parameters$calculate.empirical.p_values == TRUE){
  ## Find file with parametric bootstrap results
  all_output <- list.files(output_directory, recursive = TRUE)
  all_csv <- grep("csv", all_output, value = T)
  bs_csv_file <- paste0(output_directory, grep("empirical", grep("collated", grep("treelikeness_metrics", all_csv, value = T), value = T), value = T))
  # Open csv file
  bs_df <- read.csv(bs_csv_file, stringsAsFactors = F)
  # Add identification columns for faceting and filtering
  bs_df$rep_type <- factor(bs_df$unique_id,
                           levels = bs_df$unique_id,
                           labels = rep(c(rep("Bootstrap replicate", 100), "Alignment"), 2),
                           ordered = TRUE)
  # Add a column specifying which alignment the row is
  bs_df$alignment_id <- factor(bs_df$unique_id,
                               levels = bs_df$unique_id,
                               labels = rep(c("WEA17", "WEA17F"), each = 101),
                               ordered = TRUE)
  
  ## Calculate p-values for "WEA17" alignment
  wea17_df <- bs_df[grep("WEA17F", bs_df$unique_id, invert = T), ]
  wea17_df$unique_id[grep("rep", wea17_df$unique_id, invert = T)] <- "alignment"
  wea17_p_value_df <- calculate.all.p_values(output_df = wea17_df, 
                                             test_statistic_names = c("LM_proportion_resolved_quartets", "sCF_mean", "mean_delta_plot_value", "Cunningham_test", "tree_proportion"))
  wea17_p_value_df$dataset <- "WEA17"
  wea17_p_value_df$statistic_type  <- "p-value"
  wea17_p_value_df <- wea17_p_value_df[c("dataset", "test_statistic", "test_statistic_value", "statistic_type", "p_value")]
  names(wea17_p_value_df) <- c("dataset", "test_statistic", "test_statistic_value", "statistic_type", "statistic_value")
  # Extract Network Treelikeness Test values
  wea17_ntlt <- bs_df[which(bs_df$rep_type == "Alignment" & bs_df$alignment_id == "WEA17"),]$NetworkTreelikenessTest
  wea17_ntlt_stat_value <- length(which(bs_df$rep_type == "Bootstrap replicate" & bs_df$alignment_id == "WEA17" & bs_df$NetworkTreelikenessTest == "Treelike"))/
    length(which(bs_df$rep_type == "Bootstrap replicate" & bs_df$alignment_id == "WEA17"))
  wea17_ntlt_vector <- c("dataset" = "WEA17", "test_statistic" = "NetworkTreelikenessTest", 
                         "test_statistic_value" = wea17_ntlt, "statistic_type" = "proportion_treelike_alignments", 
                         "statistic_value" = wea17_ntlt_stat_value) 
  # Bind the tree proportion results
  wea17_p_value_df <- rbind(wea17_p_value_df, wea17_ntlt_vector)
  
  
  
  ## Calculate p-values for "WEA17F" alignment
  wea17f_df <- bs_df[grep("WEA17F", bs_df$unique_id), ]
  wea17f_df$unique_id[grep("rep", wea17f_df$unique_id, invert = T)] <- "alignment"
  wea17f_p_value_df <- calculate.all.p_values(output_df = wea17f_df, 
                                              test_statistic_names = c("LM_proportion_resolved_quartets", "sCF_mean", "mean_delta_plot_value", "Cunningham_test", "tree_proportion"))
  wea17f_p_value_df$dataset <- "WEA17F"
  wea17f_p_value_df$statistic_type  <- "p-value"
  wea17f_p_value_df <- wea17f_p_value_df[c("dataset", "test_statistic", "test_statistic_value", "statistic_type", "p_value")]
  names(wea17f_p_value_df) <- c("dataset", "test_statistic", "test_statistic_value", "statistic_type", "statistic_value")
  # Extract Network Treelikeness Test values
  wea17f_ntlt <- bs_df[which(bs_df$rep_type == "Alignment" & bs_df$alignment_id == "WEA17F"),]$NetworkTreelikenessTest
  wea17f_ntlt_stat_value <- length(which(bs_df$rep_type == "Bootstrap replicate" & bs_df$alignment_id == "WEA17F" & bs_df$NetworkTreelikenessTest == "Treelike"))/
    length(which(bs_df$rep_type == "Bootstrap replicate" & bs_df$alignment_id == "WEA17F"))
  wea17f_ntlt_vector <- c("dataset" = "WEA17F", "test_statistic" = "NetworkTreelikenessTest", 
                          "test_statistic_value" = wea17f_ntlt, "statistic_type" = "proportion_treelike_alignments", 
                          "statistic_value" = wea17f_ntlt_stat_value) 
  # Bind the tree proportion results
  wea17f_p_value_df <- rbind(wea17f_p_value_df, wea17f_ntlt_vector)
  
  ## Combine and save dataframes
  collated_p_value_df <- rbind(wea17_p_value_df, wea17f_p_value_df)
  collated_p_value_file <- paste0(output_directory, "empirical_collated_p_values.csv")
  write.csv(collated_p_value_df, file = collated_p_value_file, row.names = F)
}



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



#### 7. Apply tests for treelikeness to each gene ####
if (control_parameters$per.gene.analysis == TRUE){
  ## Split each alignment into individual genes
  # List all files in the empirical data directory
  data_files <- paste0(empirical_data_directory, list.files(empirical_data_directory, recursive = TRUE))
  data_files <- grep("genes", data_files, value = T)
  # Split Whelan 2017 alignment into separate genes
  wea17_alignment_file <- grep("alignment", grep("filtered", grep("fa", grep("WEA17.", data_files, value = T), value = T), value = T, invert = TRUE), value = T, invert = TRUE)
  wea17_partition_file <- grep("filtered", grep("partitions", grep("Whelan2017", data_files, value = T), value = T), value = T, invert = T)
  wea17_gene_directory <- paste0(dirname(wea17_alignment_file), "/")
  split.partitions(alignment_file = wea17_alignment_file, partition_file = wea17_partition_file, gene_output_directory = wea17_gene_directory)
  # Split Whelan 2017 filtered alignment into separate genes
  wea17f_alignment_file <- grep("fa", grep("WEA17F", data_files, value = T), value = T)
  wea17f_partition_file <- grep("partitions", grep("Whelan2017_filtered", data_files, value = T), value = T)
  wea17f_gene_directory <- paste0(dirname(wea17f_alignment_file), "/")
  split.partitions(alignment_file = wea17f_alignment_file, partition_file = wea17f_partition_file, gene_output_directory = wea17f_gene_directory)
  
  ## For each gene: estimate tree, then apply metrics: tree proportion, sCF, and delta plots
  
  ## For WEA17
  # Extract the list of alignments using the id (WEA17F)
  wea17_files <- paste0(replicate_alignment_directory, grep("WEA17_", all_files, value = T))
  wea17_al <- wea17_files[grep("bs_rep", basename(wea17_files), ignore.case = T, invert = T)]
  wea17_replicates <- wea17_files[grep("bs_rep", basename(wea17_files), ignore.case = T)]
  # Copy original alignment and prepare for treelikeness metric runs
  wea17_dir <- paste0(output_directory, "WEA17/")
  if (dir.exists(wea17_dir) == FALSE){ dir.create(wea17_dir) }
  wea17_copy <- paste0(wea17_dir, basename(wea17_al))
  file.copy(from = wea17_al, to = wea17_copy)
  # Calculate the tree proportion for the original alignments
  treelikeness.metrics.empirical(wea17_copy,
                                 splitstree_path = splitstree_path, 
                                 iqtree2_path = iqtree2_path, 
                                 iqtree_model = "'Q.insect+R8{0.1639,0.0315,0.1812,0.1846,0.1454,0.4618,0.1171,0.7116,0.1742,1.1628,0.1403,2.0884,0.0676,3.6840,0.0103,6.4538}'", 
                                 num_iqtree2_threads = num_cores, 
                                 sequence_format = "AA", 
                                 redo = FALSE)
  # Apply the wrapper function to calculate treelikeness of bootstrap replicates
  wea17_df <- bootstrap.wrapper(wea17_replicates, 
                                output_directory = output_directory, 
                                splitstree_path = splitstree_path, 
                                iqtree2_path = iqtree2_path, 
                                iqtree_model = "'Q.insect+R8{0.1639,0.0315,0.1812,0.1846,0.1454,0.4618,0.1171,0.7116,0.1742,1.1628,0.1403,2.0884,0.0676,3.6840,0.0103,6.4538}'", 
                                num_iqtree2_threads = "10", 
                                sequence_format = "AA", 
                                redo = FALSE, 
                                number_parallel_cores = mclapply_num_cores)
  
  ## For WEA17_filtered
  # Extract the list of alignments using the id (WEA17F)
  wea17f_files <- paste0(replicate_alignment_directory, grep("WEA17F_", all_files, value = T))
  wea17f_al <- wea17f_files[grep("bs_rep", basename(wea17f_files), ignore.case = T, invert = T)]
  wea17f_replicates <- wea17f_files[grep("bs_rep", basename(wea17f_files), ignore.case = T)]
  # Copy original alignment and prepare for treelikeness metric runs
  wea17f_dir <- paste0(output_directory, "WEA17F/")
  if (dir.exists(wea17f_dir) == FALSE){ dir.create(wea17f_dir) }
  wea17f_copy <- paste0(wea17f_dir, basename(wea17f_al))
  file.copy(from = wea17f_al, to = wea17f_copy)
  # Calculate the tree proportion for the original alignments
  treelikeness.metrics.empirical(wea17f_copy,
                                 splitstree_path = splitstree_path, 
                                 iqtree2_path = iqtree2_path, 
                                 iqtree_model = "'Q.insect+I{0.0659}+R6{0.2004,0.1133,0.1734,0.3827,0.1439,0.6654,0.2108,1.1677,0.1549,2.2361,0.0508,4.3870}'", 
                                 num_iqtree2_threads = num_cores, 
                                 sequence_format = "AA", 
                                 redo = FALSE)
  
  # Apply the wrapper function
  wea17f_df <- bootstrap.wrapper(wea17f_replicates, 
                                 output_directory = output_directory, 
                                 splitstree_path = splitstree_path, 
                                 iqtree2_path = iqtree2_path, 
                                 iqtree_model = "'Q.insect+I{0.0659}+R6{0.2004,0.1133,0.1734,0.3827,0.1439,0.6654,0.2108,1.1677,0.1549,2.2361,0.0508,4.3870}'", 
                                 num_iqtree2_threads = "10", 
                                 sequence_format = "AA", 
                                 redo = FALSE, 
                                 number_parallel_cores = mclapply_num_cores)
}