## Go back to how you are scaling the trees and see how it's possible that the
#   tree depths can be so much larger than the target.



#### 1. Prepare packages ####
library(treestats)
library(ape)
library(phytools)
library(parallel)
library(ggplot2)
library(tidyr)
library(dplyr)
library(forcats)



#### 2. Prepare functions ####
check.tree.scaling <- function(row_number, params_df) {
  ## Generate random trees and check scaling
  print(row_number)
  row <- params_df[row_number, ]
  rownames(row) <- NULL
  ## Test 1: random trees
  # NOTE:  Code taken from function "generate.random,trees" in 'code/func_simulating_alignments.R'
  # Generate random trees (these are not ultrametric!)
  random_trees <- rmtree(row$num_trees, row$num_taxa)
  # Scale the depth of each tree (so the total depth is set according to the tree_depth parameter)
  scaled_trees <- random_trees
  for (i in 1:length(scaled_trees)) {
    t <- scaled_trees[[i]]
    t$edge.length <- t$edge.length * (row$tree_depth / max(branching.times(t)))
    scaled_trees[[i]] <- t
  }
  # Calculate stats for scaled random trees
  st_stats <- apply.treestats(scaled_trees)
  names(st_stats) <- paste0("st_rt_", names(st_stats))
  ## Test 2: random (scaled) trees
  # Generate random trees (these are not ultrametric!)
  random_ext_trees <- rmtree(row$num_trees, row$num_taxa)
  # Extend tips so each tree is ultrametric
  scaled_ext_trees <- lapply(random_ext_trees, function(t) {
    force.ultrametric(t, method = "extend", message = FALSE)
  })
  # Scale the depth of each tree (so the total depth is set according to the tree_depth parameter)
  for (i in 1:length(scaled_ext_trees)) {
    t <- scaled_ext_trees[[i]]
    t$edge.length <- t$edge.length * (row$tree_depth / max(branching.times(t)))
    scaled_ext_trees[[i]] <- t
  }
  # Calculate stats for scaled random trees
  st_ext_stats <- apply.treestats(scaled_ext_trees)
  names(st_ext_stats) <- paste0("st_ext_", names(st_ext_stats))
  ## Test 3: coalescent random trees
  # Generate random trees (these are coalescent and ultrametric)
  random_coal_trees <- lapply(1:row$num_trees, function(i) {
    rcoal(n = row$num_taxa)
  })
  # Scale the depth of each tree (so the total depth is set according to the tree_depth parameter)
  scaled_coal_trees <- random_coal_trees
  for (i in 1:length(scaled_coal_trees)) {
    t <- scaled_coal_trees[[i]]
    t$edge.length <- t$edge.length * (row$tree_depth / max(branching.times(t)))
    scaled_coal_trees[[i]] <- t
  }
  # Calculate stats for scaled random trees
  coal_st_stats <- apply.treestats(scaled_coal_trees)
  names(coal_st_stats) <- paste0("st_coal_", names(coal_st_stats))
  ## Combine results
  op <- c(st_stats, st_ext_stats, coal_st_stats)
  op_df <- as.data.frame(matrix(data = op, nrow = 1, ncol = length(op), byrow = TRUE))
  names(op_df) <- names(op)
  op_df <- cbind(row, op_df)
  return(op_df)
}

apply.treestats <- function(multiphylo) {
  max_branching_times <- unlist(lapply(multiphylo, function(t) {
    max(treestats::branching_times(t))
  }))
  # Tree height = the maximum branching time plus the root branch length
  tree_height <- unlist(lapply(multiphylo, function(t) {
    treestats::tree_height(t)
  }))
  # Crown age = the maximum branching time
  crown_age <- unlist(lapply(multiphylo, function(t) {
    treestats::crown_age(t)
  }))
  stats <- c(summary(max_branching_times),
             summary(tree_height),
             summary(crown_age))
  names(stats) <- paste0(rep(c(
    "max_branching_times", "tree_height", "crown_age"
  ), each = 6), "_", rep(c(
    "min", "1st_qu", "median", "mean", "3rd_qu", "max"
  ), 3))
  return(stats)
}

divisors <- function(x){
  # Function to find all divisors of x using modulo division
  y <- 1:x
  # If remainder is 0, then y is a divisor of x. Keep all divisors
  d = y[x%%y == 0]
  return(d)
}



#### 3. Create dataframe for testing tree scaling ####
ts_params <- expand.grid(
  "num_reps" = 1:20,
  "num_taxa" = c(5, 10, 20, 50, 100),
  "num_trees" = divisors(10000),
  "tree_depth" = c(0.01, 0.1, 1),
  "alisim_gene_models" = "JC",
  "alisim_gene_tree_length" = NA,
  "total_alignment_length" = 10000,
  "sequence_type" = "DNA"
)
# Add a unique identifier (uid) of the form:
#   ts_{num. trees}_{num. taxa}_{rep. num.}_{tree_depth}
ts_params$uid <- paste0(
  "ts_",
  sprintf("%05d", ts_params$num_trees),
  "_",
  sprintf("%04d", ts_params$num_taxa),
  "_",
  sprintf("%03d", ts_params$num_reps),
  "_",
  ts_params$tree_depth
)



#### 4. Generate random trees and return scaled tree depth ####
ts_stats_file <- "supp_output/check_01_tree_scaling.csv"
if (file.exists(ts_stats_file) == FALSE) {
  ts_list <- mclapply(
    1:nrow(ts_params),
    check.tree.scaling,
    params_df = ts_params,
    mc.cores = 10
  )
  ts_stats <- as.data.frame(do.call(rbind, ts_list))
  write.csv(ts_stats, file = ts_stats_file, row.names = FALSE)
} else {
  ts_stats <- read.csv(ts_stats_file)
}



#### 5. Plot tree depth for different random tree methods ####
## Key:
# st_rt = scaled trees, random trees (generated with ape::rmtree)
# st_ext = scaled trees, random trees (generated with ape::rmtree),
#          made ultrametric (with phytools::force.ultrametric, method = "extend")
# st_coal = scaled trees, coalescent trees (generated with ape::rcoal)

# Plot 1: Plot histogram of median maximum branching times
p_df <- ts_stats %>%
  select(
    num_reps:uid,
    st_rt_max_branching_times_median,
    st_ext_max_branching_times_median,
    st_coal_max_branching_times_median
  ) %>%
  pivot_longer(cols = starts_with("st_"),
               names_to = "variable",
               values_to = "value") %>%
  mutate(tree_method = unlist(lapply(strsplit(variable, "_"), function(x) {
    paste(x[1:2], collapse = "_")
  })),
  tree_stat = unlist(lapply(strsplit(variable, "_"), function(x) {
    paste(x[3:length(x)], collapse = "_")
  }))) %>%
  mutate(
    tree_method = case_when(
      tree_method == "st_rt" ~ "rtree",
      tree_method == "st_ext" ~ "rtree_extended",
      tree_method == "st_coal" ~ "rcoal"
    )
  ) %>%
  mutate(
    tree_method = as.factor(tree_method) %>%
      fct_relevel(., c("rtree", "rtree_extended", "rcoal")) %>%
      fct_recode(
        .,
        "Random tree" = "rtree",
        "Random tree\n(extended tips)" = "rtree_extended",
        "Random coalescent\ntree" = "rcoal"
      )
  )
p <- ggplot(p_df, aes(x = value)) +
  facet_grid(tree_method ~ tree_depth, scales = "free") +
  geom_histogram() +
  xlab("Median maximum branching time") +
  ylab("Count") +
  labs(title = "Median maximum branching time for all random trees\nin a single simulation replicate")
ggsave(filename = "supp_plots/tree_scaling_check_hist_max_branching_time.pdf",
       plot = p)
ggsave(filename = "supp_plots/tree_scaling_check_hist_max_branching_time.png",
       plot = p)

# Plot 2: Plot histogram of median tree height
p_df <- ts_stats %>%
  select(
    num_reps:uid,
    st_rt_tree_height_median,
    st_ext_tree_height_median,
    st_coal_tree_height_median
  ) %>%
  pivot_longer(cols = starts_with("st_"),
               names_to = "variable",
               values_to = "value") %>%
  mutate(tree_method = unlist(lapply(strsplit(variable, "_"), function(x) {
    paste(x[1:2], collapse = "_")
  })),
  tree_stat = unlist(lapply(strsplit(variable, "_"), function(x) {
    paste(x[3:length(x)], collapse = "_")
  }))) %>%
  mutate(
    tree_method = case_when(
      tree_method == "st_rt" ~ "rtree",
      tree_method == "st_ext" ~ "rtree_extended",
      tree_method == "st_coal" ~ "rcoal"
    )
  ) %>%
  mutate(
    tree_method = as.factor(tree_method) %>%
      fct_relevel(., c("rtree", "rtree_extended", "rcoal")) %>%
      fct_recode(
        .,
        "Random tree" = "rtree",
        "Random tree\n(extended tips)" = "rtree_extended",
        "Random coalescent\ntree" = "rcoal"
      )
  )
p <- ggplot(p_df, aes(x = value)) +
  facet_grid(tree_method ~ tree_depth, scales = "free") +
  geom_histogram() +
  xlab("Median tree height") +
  ylab("Count") +
  labs(title = "Median tree height for all random trees\nin a single simulation replicate",
       subtitle = "Tree height = max. branching time plus root branch length")
ggsave(filename = "supp_plots/tree_scaling_check_hist_tree_height.pdf",
       plot = p)
ggsave(filename = "supp_plots/tree_scaling_check_hist_tree_height.png",
       plot = p)

# Plot 3: Plot histogram of median crown age
p_df <- ts_stats %>%
  select(
    num_reps:uid,
    st_rt_crown_age_median,
    st_ext_crown_age_median,
    st_coal_crown_age_median
  ) %>%
  pivot_longer(cols = starts_with("st_"),
               names_to = "variable",
               values_to = "value") %>%
  mutate(tree_method = unlist(lapply(strsplit(variable, "_"), function(x) {
    paste(x[1:2], collapse = "_")
  })),
  tree_stat = unlist(lapply(strsplit(variable, "_"), function(x) {
    paste(x[3:length(x)], collapse = "_")
  }))) %>%
  mutate(
    tree_method = case_when(
      tree_method == "st_rt" ~ "rtree",
      tree_method == "st_ext" ~ "rtree_extended",
      tree_method == "st_coal" ~ "rcoal"
    )
  ) %>%
  mutate(
    tree_method = as.factor(tree_method) %>%
      fct_relevel(., c("rtree", "rtree_extended", "rcoal")) %>%
      fct_recode(
        .,
        "Random tree" = "rtree",
        "Random tree\n(extended tips)" = "rtree_extended",
        "Random coalescent\ntree" = "rcoal"
      )
  )
p <- ggplot(p_df, aes(x = value)) +
  facet_grid(tree_method ~ tree_depth, scales = "free") +
  geom_histogram() +
  xlab("Median crown age") +
  ylab("Count") +
  labs(title = "Median crown age for all random trees\nin a single simulation replicate",
       subtitle = "Crown age = max. branching time")
ggsave(filename = "supp_plots/tree_scaling_check_hist_crown_age.pdf",
       plot = p)
ggsave(filename = "supp_plots/tree_scaling_check_hist_crown_age.png",
       plot = p)

# Plot 4: Plot tree depth against median maximum branching times
p_df <- ts_stats %>%
  select(
    num_reps:uid,
    st_rt_max_branching_times_median,
    st_ext_max_branching_times_median,
    st_coal_max_branching_times_median
  ) %>%
  pivot_longer(cols = starts_with("st_"),
               names_to = "variable",
               values_to = "value") %>%
  mutate(tree_method = unlist(lapply(strsplit(variable, "_"), function(x) {
    paste(x[1:2], collapse = "_")
  })),
  tree_stat = unlist(lapply(strsplit(variable, "_"), function(x) {
    paste(x[3:length(x)], collapse = "_")
  }))) %>%
  mutate(
    tree_method = case_when(
      tree_method == "st_rt" ~ "rtree",
      tree_method == "st_ext" ~ "rtree_extended",
      tree_method == "st_coal" ~ "rcoal"
    )
  ) %>%
  mutate(
    tree_method = as.factor(tree_method) %>%
      fct_relevel(., c("rtree", "rtree_extended", "rcoal")) %>%
      fct_recode(
        .,
        "Random tree" = "rtree",
        "Random tree\n(extended tips)" = "rtree_extended",
        "Random coalescent\ntree" = "rcoal"
      )
  ) %>%
  mutate(
    tree_depth = as.factor(as.character(tree_depth)) %>%
      fct_relevel(., c("0.01", "0.1", "1")) %>%
      fct_recode(., "0.01" = "0.01", "0.10" = "0.1", "1.00" = "1")
  )
p <- ggplot(p_df, aes(x = tree_method, y = value, colour = tree_depth)) +
  geom_boxplot() +
  xlab("Tree simulation method") +
  ylab("Median maximum branching time") +
  labs(title = "Median maximum branching time for all random trees\ngenerated by three methods") +
  scale_colour_viridis_d(option = "H") +
  guides(colour = guide_legend(title = "Tree depth")) +
  theme_bw()
ggsave(filename = "supp_plots/tree_scaling_check_boxplot_max_branching_time.pdf",
       plot = p)
ggsave(filename = "supp_plots/tree_scaling_check_boxplot_max_branching_time.png",
       plot = p)

# Plot 5: Plot tree depth against median tree height
p_df <- ts_stats %>%
  select(
    num_reps:uid,
    st_rt_tree_height_median,
    st_ext_tree_height_median,
    st_coal_tree_height_median
  ) %>%
  pivot_longer(cols = starts_with("st_"),
               names_to = "variable",
               values_to = "value") %>%
  mutate(tree_method = unlist(lapply(strsplit(variable, "_"), function(x) {
    paste(x[1:2], collapse = "_")
  })),
  tree_stat = unlist(lapply(strsplit(variable, "_"), function(x) {
    paste(x[3:length(x)], collapse = "_")
  }))) %>%
  mutate(
    tree_method = case_when(
      tree_method == "st_rt" ~ "rtree",
      tree_method == "st_ext" ~ "rtree_extended",
      tree_method == "st_coal" ~ "rcoal"
    )
  ) %>%
  mutate(
    tree_method = as.factor(tree_method) %>%
      fct_relevel(., c("rtree", "rtree_extended", "rcoal")) %>%
      fct_recode(
        .,
        "Random tree" = "rtree",
        "Random tree\n(extended tips)" = "rtree_extended",
        "Random coalescent\ntree" = "rcoal"
      )
  ) %>%
  mutate(
    tree_depth = as.factor(as.character(tree_depth)) %>%
      fct_relevel(., c("0.01", "0.1", "1")) %>%
      fct_recode(., "0.01" = "0.01", "0.10" = "0.1", "1.00" = "1")
  )
p <- ggplot(p_df, aes(x = tree_method, y = value, colour = tree_depth)) +
  geom_boxplot() +
  xlab("Tree simulation method") +
  ylab("Median tree height") +
  labs(title = "Median tree height for all random trees generated by three methods",
       subtitle = "Tree height = max. branching time plus root branch length") +
  scale_colour_viridis_d(option = "H") +
  guides(colour = guide_legend(title = "Tree depth")) +
  theme_bw()
ggsave(filename = "supp_plots/tree_scaling_check_boxplot_tree_height.pdf",
       plot = p)
ggsave(filename = "supp_plots/tree_scaling_check_boxplot_tree_height.png",
       plot = p)

# Plot 6: Plot tree depth against median crown age
p_df <- ts_stats %>%
  select(
    num_reps:uid,
    st_rt_crown_age_median,
    st_ext_crown_age_median,
    st_coal_crown_age_median
  ) %>%
  pivot_longer(cols = starts_with("st_"),
               names_to = "variable",
               values_to = "value") %>%
  mutate(tree_method = unlist(lapply(strsplit(variable, "_"), function(x) {
    paste(x[1:2], collapse = "_")
  })),
  tree_stat = unlist(lapply(strsplit(variable, "_"), function(x) {
    paste(x[3:length(x)], collapse = "_")
  }))) %>%
  mutate(
    tree_method = case_when(
      tree_method == "st_rt" ~ "rtree",
      tree_method == "st_ext" ~ "rtree_extended",
      tree_method == "st_coal" ~ "rcoal"
    )
  ) %>%
  mutate(
    tree_method = as.factor(tree_method) %>%
      fct_relevel(., c("rtree", "rtree_extended", "rcoal")) %>%
      fct_recode(
        .,
        "Random tree" = "rtree",
        "Random tree\n(extended tips)" = "rtree_extended",
        "Random coalescent\ntree" = "rcoal"
      )
  ) %>%
  mutate(
    tree_depth = as.factor(as.character(tree_depth)) %>%
      fct_relevel(., c("0.01", "0.1", "1")) %>%
      fct_recode(., "0.01" = "0.01", "0.10" = "0.1", "1.00" = "1")
  )
p <- ggplot(p_df, aes(x = tree_method, y = value, colour = tree_depth)) +
  geom_boxplot() +
  xlab("Tree simulation method") +
  ylab("Median crown age") +
  labs(title = "Median crown age for all random trees generated by three methods",
       subtitle = "Crown age = max. branching time") +
  scale_colour_viridis_d(option = "H") +
  guides(colour = guide_legend(title = "Tree depth")) +
  theme_bw()
ggsave(filename = "supp_plots/tree_scaling_check_boxplot_crown_age.pdf",
       plot = p)
ggsave(filename = "supp_plots/tree_scaling_check_boxplot_crown_age.png",
       plot = p)

# Plot 7: Plot number of taxa against median crown age
p_df <- ts_stats %>%
  select(
    num_reps:uid,
    st_rt_crown_age_median
  ) %>%
  pivot_longer(cols = starts_with("st_"),
               names_to = "variable",
               values_to = "value") %>%
  mutate(tree_method = unlist(lapply(strsplit(variable, "_"), function(x) {
    paste(x[1:2], collapse = "_")
  })),
  tree_stat = unlist(lapply(strsplit(variable, "_"), function(x) {
    paste(x[3:length(x)], collapse = "_")
  }))) %>%
  mutate(
    tree_method = case_when(
      tree_method == "st_rt" ~ "rtree",
      tree_method == "st_ext" ~ "rtree_extended",
      tree_method == "st_coal" ~ "rcoal"
    )
  ) %>%
  mutate(
    tree_method = as.factor(tree_method) %>%
      fct_relevel(., c("rtree")) %>%
      fct_recode(
        .,
        "Random tree" = "rtree"
      )
  ) %>%
  mutate(
    tree_depth = as.factor(as.character(tree_depth)) %>%
      fct_relevel(., c("0.01", "0.1", "1")) %>%
      fct_recode(., "0.01" = "0.01", "0.10" = "0.1", "1.00" = "1")
  ) %>%
  mutate(
    num_taxa = as.factor(as.character(num_taxa)) %>%
      fct_relevel(., c("5", "10", "20", "50", "100"))
  )
p <- ggplot(p_df, aes(x = tree_method, y = value, colour = num_taxa)) +
  facet_wrap(tree_depth ~ ., scales = "free" ) +
  geom_boxplot() +
  xlab("Tree simulation method") +
  ylab("Median crown age") +
  labs(title = "Median crown age for trees with different numbers of taxa",
       subtitle = "Crown age = max. branching time") +
  scale_colour_viridis_d(option = "H") +
  guides(colour = guide_legend(title = "Num. taxa")) +
  theme_bw()
ggsave(filename = "supp_plots/tree_scaling_check_boxplot_numTaxa_crownAge.pdf",
       plot = p)
ggsave(filename = "supp_plots/tree_scaling_check_boxplot_numTaxa_crownAge.png",
       plot = p)
