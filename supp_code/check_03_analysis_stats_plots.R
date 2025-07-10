## Plotting analysis statistics

## 1. Prepare input
# Statistics files
b1_df_file <- "supp_output/expB1_check_02_analysis_stats_collated.csv"
# Output plot directory
plot_dir <- "supp_plots"
if (dir.exists(plot_dir) == FALSE){dir.create(plot_dir)}


## 2. Prepare packages
library(ggplot2)
library(dplyr)
library(tidyr)


## 3. Reformat dataframes
b1_df_raw <- read.csv(b1_df_file)


## 4. Plot IQ-Tree statistics
# b1_df <- b1_df_raw %>%
#   select(!(
#     starts_with("random_trees") |
#       starts_with("al_") |
#       starts_with("iqtree_tree") |
#       starts_with("quartet_treeness")
#   )) %>%
#   select(!all_of(
#     c(
#       "row_id",
#       "tree_file",
#       "partition_file",
#       "output_alignment_file"
#     )
#   )) %>%
#   pivot_longer(
#     cols = num_sites:percent_parsimony_informative_sites,
#     names_to = "variable",
#     values_to = "value"
#   )

## Histogram: number of invariant sites, faceted by tree depth
p_df <- b1_df_raw %>%
  select(!(
    starts_with("random_trees") |
      starts_with("al_") |
      starts_with("iqtree_tree") |
      starts_with("quartet_treeness")
  )) %>%
  select(!all_of(
    c(
      "row_id",
      "alisim_gene_tree_length",
      "tree_file",
      "partition_file",
      "output_alignment_file"
    )
  ))
p <- ggplot(p_df, aes(x = number_constant_sites)) +
  facet_wrap(vars(tree_depth)) +
  geom_histogram(binwidth = 500) +
  xlab("Number of constant sites") +
  ylab("Count") +
  labs(title = "B2: Constant sites for three tree depths")
ggsave(
  filename = paste0(plot_dir, "/B2_hist_constant_sites.pdf"),
  plot = p,
  device = "pdf"
)

## Histogram: number of parsimoney-informative sites, faceted by tree depth
p_df <- b1_df_raw %>%
  select(!(
    starts_with("random_trees") |
      starts_with("al_") |
      starts_with("iqtree_tree") |
      starts_with("quartet_treeness")
  )) %>%
  select(!all_of(
    c(
      "row_id",
      "alisim_gene_tree_length",
      "tree_file",
      "partition_file",
      "output_alignment_file"
    )
  ))
p <- ggplot(p_df, aes(x = number_parsimony_informative_sites)) +
  facet_wrap(vars(tree_depth)) +
  geom_histogram(binwidth = 500) +
  xlab("Number of parsimony-informative sites") +
  ylab("Count") +
  labs(title = "B2: Parsimony-informative sites for three tree depths")
ggsave(
  filename = paste0(plot_dir, "/B2_hist_parsimony-informative_sites.pdf"),
  plot = p,
  device = "pdf"
)


## Histogram: raw pairwise distances for alignment, faceted by tree depth
p_df <- b1_df_raw %>%
  select(!(
    starts_with("random_trees") |
      starts_with("iqtree_tree") |
      starts_with("quartet_treeness") |
      starts_with("number") |
      starts_with("percent") |
      starts_with("al_jc")
  )) %>%
  select (! ends_with("quartets")) %>%
  select(!all_of(
    c(
      "row_id",
      "alisim_gene_tree_length",
      "tree_file",
      "partition_file",
      "output_alignment_file"
    )
  ))
# Median pwd
p <- ggplot(p_df, aes(x = al_raw_pairwise_distance_median)) +
  facet_wrap(vars(tree_depth)) +
  geom_histogram(bins = 25) +
  xlab("Median raw pairwise distance (alignment)") +
  ylab("Count") +
  labs(title = "B2: Alignment pairwise distances (model = raw) for three tree depths")
ggsave(
  filename = paste0(plot_dir, "/B2_hist_median_alnmt_raw_pwd.pdf"),
  plot = p,
  device = "pdf"
)
# Min pwd
p <- ggplot(p_df, aes(x = al_raw_pairwise_distance_min)) +
  facet_wrap(vars(tree_depth)) +
  geom_histogram(bins = 25) +
  xlab("Min. raw pairwise distance (alignment)") +
  ylab("Count") +
  labs(title = "B2: Alignment pairwise distances (model = raw) for three tree depths")
ggsave(
  filename = paste0(plot_dir, "/B2_hist_min_alnmt_raw_pwd.pdf"),
  plot = p,
  device = "pdf"
)
# Max pwd
p <- ggplot(p_df, aes(x = al_raw_pairwise_distance_max)) +
  facet_wrap(vars(tree_depth)) +
  geom_histogram(bins = 25) +
  xlab("Max. raw pairwise distance (alignment)") +
  ylab("Count") +
  labs(title = "B2: Alignment pairwise distances (model = raw) for three tree depths")
ggsave(
  filename = paste0(plot_dir, "/B2_hist_max_alnmt_raw_pwd.pdf"),
  plot = p,
  device = "pdf"
)

## Scatterplot: median random tree max branching time vs iqtree tree max branching time
p_df <- b1_df_raw %>%
  select(all_of(
    c(
      "uid",
      "num_reps",
      "num_taxa",
      "num_trees",
      "tree_depth",
      "alisim_gene_models",
      "total_alignment_length",
      "random_trees_max_branching_time_median",
      "iqtree_tree_max_branching_time"
    )
  ))
# Maximum branching time for random trees
p <- ggplot(p_df, aes(x = random_trees_max_branching_time_median)) +
  facet_wrap(vars(tree_depth), scales = "free") +
  geom_histogram() +
  xlab("Median maximum branching time across all random trees") +
  ylab("Count") +
  labs(title = "B2: Randomly generated trees - maximum branching time for three tree depths")
ggsave(
  filename = paste0(plot_dir, "/B2_hist_random_tree_max_branching_time.pdf"),
  plot = p,
  device = "pdf"
)
# Maximum branching time for IQ-Tree tree
p <- ggplot(p_df, aes(x = iqtree_tree_max_branching_time)) +
  facet_wrap(vars(tree_depth), scales = "free") +
  geom_histogram() +
  xlab("Maximum branching time") +
  ylab("Count") +
  labs(title = "B2: IQ-Tree ML Tree - maximum branching time for three tree depths")
ggsave(
  filename = paste0(plot_dir, "/B2_hist_iqtree_tree_max_branching_time.pdf"),
  plot = p,
  device = "pdf"
)
# Plot branching times against each other (random tree vs iqtree)
p <- ggplot(p_df, aes(x = random_trees_max_branching_time_median, y = iqtree_tree_max_branching_time)) +
  facet_wrap(vars(tree_depth), scales = "free") +
  geom_point() +
  xlab("Median maximum branching time across all random trees") +
  ylab("Maximum branching time for inferred IQ-Tree tree") +
  labs(title = "B2: comparing maximum branching time for random trees and IQ-tree tree")
ggsave(
  filename = paste0(plot_dir, "/B2_scatter_compare_max_branching_times.pdf"),
  plot = p,
  device = "pdf"
)
