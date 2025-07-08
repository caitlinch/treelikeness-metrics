library(ape)
library(treestats)

check.expA1.analysis.stats <- function(simulation_directory){
  ## Calculate statistics on the alignment and trees for each replicate in expA1
  # List all files in this directory
  rep_files <- list.files(simulation_directory)
  # Open simulation parameters
  rep_params <- read.csv(paste0(
    simulation_directory,
    grep("parameters.csv", rep_files, value = TRUE)
  ))
  # Open NTLT results
  rep_results <- read.csv(paste0(
    simulation_directory,
    grep("_NTLT_results.csv", rep_files, value = TRUE)
  ))
  rep_results <- rep_results[, which(names(rep_results) == "identical_sequences_present"):ncol(rep_results)]
  # Extract statistics about the alignment
  rep_al_file <- paste0(simulation_directory,
                        grep(
                          "output_alignment.fa.",
                          grep("output_alignment.fa", rep_files, value = TRUE),
                          invert = TRUE,
                          value = TRUE
                        ))
  rep_pwd <- check.pairwise.distances(rep_al_file)
  # Extract tree statistics for simulated trees
  rep_trees_file <- paste0(
    simulation_directory,
    grep("random_trees.phy", rep_files, value = TRUE)
  )
  rep_trees <- ape::read.tree(rep_trees_file)
  rep_trees_values <- check.tree.stats(rep_trees)
  # Extract information about IQ-Tree tree
  rep_iqtree_treefile <- paste0(
    simulation_directory,
    grep("output_alignment.fa.treefile", rep_files, value = TRUE)
  )
  rep_iqtree_tree_stats <- check.tree.stats(
    read.tree(rep_iqtree_treefile)
  )
  names(rep_iqtree_tree_stats) <- paste0(
    "iqtree_tree_",
    names(rep_iqtree_tree_stats)
  )
  # Extract information from the IQ-Tree .iqtree file
  rep_iqtree_file <- paste0(
    simulation_directory,
    grep("output_alignment.fa.iqtree", rep_files, value = TRUE)
  )
  rep_iqtree_values <- check.iqtree.file(rep_iqtree_file)
  # Collate function output
  op_data <- c(rep_pwd,
               rep_trees_values,
               rep_iqtree_values,
               rep_iqtree_tree_stats)
  op_df <- as.data.frame(matrix(
    data = op_data,
    nrow = 1,
    ncol = length(op_data),
    byrow = TRUE
  ))
  names(op_df) <- c(
    names(rep_pwd),
    names(rep_trees_values),
    names(rep_iqtree_values),
    names(rep_iqtree_tree_stats)
  )
  # Collate with rep_params object
  op_df <- cbind(rep_params, rep_results, op_df)
  # Save dataframe
  write.csv(
    op_df,
    file = paste0(
      simulation_directory,
      op_df$uid,
      ".check_02_analysis_stats.csv"),
    row.names = FALSE
  )
}


check.expA2.analysis.stats <- function(simulation_directory){
  ## Calculate statistics on the alignment and trees for each replicate in expA1
  # List all files in this directory
  rep_files <- list.files(simulation_directory)
  if (TRUE %in% grepl("_NTLT_results.csv", rep_files)){
    NTLT_file = paste0(simulation_directory, grep("_NTLT_results.csv", rep_files, value = TRUE))
    if (file.info(NTLT_file)[["size"]] != 0){
      # Open simulation parameters
      rep_params <- read.csv(paste0(
        simulation_directory,
        grep("parameters.csv", rep_files, value = TRUE)
      ))
      # Open NTLT results
      rep_results <- read.csv(paste0(
        simulation_directory,
        grep("_NTLT_results.csv", rep_files, value = TRUE)
      ))
      rep_results <- rep_results[, which(names(rep_results) == "identical_sequences_present"):ncol(rep_results)]
      # Extract statistics about the alignment
      rep_al_file <- paste0(simulation_directory,
                            grep(
                              "output_alignment.fa.",
                              grep("output_alignment.fa", rep_files, value = TRUE),
                              invert = TRUE,
                              value = TRUE
                            ))
      rep_pwd <- check.pairwise.distances(rep_al_file)
      # Extract tree statistics for simulated trees
      rep_trees_file <- paste0(simulation_directory, grep("random_trees.phy", rep_files, value = TRUE))
      rep_trees <- ape::read.tree(rep_trees_file)
      if (is(rep_trees, "phylo")){
        check_trees_values <- check.tree.stats(rep_trees)
        rep_trees_values <- c(
          check_trees_values,
          rep(NA, 144)
        )
        names(rep_trees_values) <-     paste0(
          rep(paste0("tree_", 1:10, "_"), each = 16),
          rep(names(check_trees_values), 10)
        )
      } else if (is(rep_trees, "multiPhylo")){
        check_trees_values <- lapply(rep_trees, check.tree.stats)
        present_trees <- length(check_trees_values)
        missing_trees <- (present_trees + 1):10
        for (i in missing_trees){
          check_trees_values[[i]] <- c(
            "mean_nearest_taxon_distance" = NA,
            "mean_pairwise_distance" = NA,
            "mean_branch_length" = NA,
            "var_branch_length" = NA,
            "mean_branch_length_external" = NA,
            "var_branch_length_external" = NA,
            "mean_branch_length_internal" = NA,
            "var_branch_length_internal" = NA,
            "average_vertex_depth" = NA,
            "max_width" = NA,
            "max_depth" = NA,
            "max_ladder" = NA,
            "max_branching_time" = NA,
            "crown_age" = NA,
            "tree_height" = NA,
            "treeness" = NA
          )
        }
        for (i in 1:length(check_trees_values)){
          names(check_trees_values[[i]]) <- paste0("tree_", i, "_",  names(check_trees_values[[i]]))
        }
        rep_trees_values <- unlist(check_trees_values)
      }
      # need output = rep_trees_values, with 10 trees worth of data
      # Extract information about IQ-Tree tree
      rep_iqtree_treefile <- paste0(
        simulation_directory,
        grep("output_alignment.fa.treefile", rep_files, value = TRUE)
      )
      rep_iqtree_tree_stats <- check.tree.stats(
        read.tree(rep_iqtree_treefile)
      )
      names(rep_iqtree_tree_stats) <- paste0("iqtree_tree_", names(rep_iqtree_tree_stats))
      # Extract information from the IQ-Tree .iqtree file
      rep_iqtree_file <- paste0(
        simulation_directory,
        grep("output_alignment.fa.iqtree", rep_files, value = TRUE)
      )
      rep_iqtree_values <- check.iqtree.file(rep_iqtree_file)
      # Collate function output
      op_data <- c(rep_pwd,
                   rep_trees_values,
                   rep_iqtree_values,
                   rep_iqtree_tree_stats)
      op_df <- as.data.frame(matrix(
        data = op_data,
        nrow = 1,
        ncol = length(op_data),
        byrow = TRUE
      ))
      names(op_df) <- c(
        names(rep_pwd),
        names(rep_trees_values),
        names(rep_iqtree_values),
        names(rep_iqtree_tree_stats)
      )
      # Collate with rep_params object
      op_df <- cbind(rep_params, rep_results, op_df)
      # Save dataframe
      write.csv(
        op_df,
        file = paste0(simulation_directory, op_df$uid, ".check_02_analysis_stats.csv"),
        row.names = FALSE
      )
    }
  }
}


check.pairwise.distances <- function(alignment_file){
  # Check pairwise distances in alignment
  alnmt <- read.FASTA(alignment_file)
  jc_matrix <- dist.dna(alnmt, model = "JC")
  jc_vals <- summary(as.numeric(jc_matrix))
  raw_matrix <- dist.dna(alnmt, model = "raw")
  raw_vals <- summary(as.numeric(raw_matrix))
  al_stats <- c(
    "al_raw_pairwise_distance_min" = raw_vals[["Min."]],
    "al_raw_pairwise_distance_1st_qu" = raw_vals[["1st Qu."]],
    "al_raw_pairwise_distance_median" = raw_vals[["Median"]],
    "al_raw_pairwise_distance_mean" = raw_vals[["Mean"]],
    "al_raw_pairwise_distance_3rd_qu" = raw_vals[["3rd Qu."]],
    "al_raw_pairwise_distance_max" = raw_vals[["Max."]],
    "al_raw_pairwise_distance_sd" = sd(as.numeric(raw_matrix)),
    "al_jc_pairwise_distance_min" = jc_vals[["Min."]],
    "al_jc_pairwise_distance_1st_qu" = jc_vals[["1st Qu."]],
    "al_jc_pairwise_distance_median" = jc_vals[["Median"]],
    "al_jc_pairwise_distance_mean" = jc_vals[["Mean"]],
    "al_jc_pairwise_distance_3rd_qu" = jc_vals[["3rd Qu."]],
    "al_jc_pairwise_distance_max" = jc_vals[["Max."]],
    "al_jc_pairwise_distance_sd" = sd(as.numeric(jc_matrix))
  )
  return(al_stats)
}


check.iqtree.file <- function(iqtree_file){
  # Extract details about alignment sites from IQ-Tree .iqtree file
  iqtree_lines <- readLines(iqtree_file)
  iqtree_check_1 <- iqtree_lines[grep("Input data:", iqtree_lines)]
  iqtree_check_2 <- iqtree_lines[grep("Number of constant sites:", iqtree_lines)]
  iqtree_check_3 <- iqtree_lines[grep("Number of invariant \\(constant or ambiguous constant\\) sites:",
                                      iqtree_lines)]
  iqtree_check_4 <- iqtree_lines[grep("Number of parsimony informative sites:", iqtree_lines)]
  iqtree_check_5 <- iqtree_lines[grep("Number of distinct site patterns:", iqtree_lines)]
  iqtree_values <- c(
    "num_sites" = as.numeric(gsub(
      " ", "", strsplit(strsplit(iqtree_check_1, "with")[[1]][2], "nucleotide")[[1]][1]
    )),
    "number_constant_sites" = as.numeric(gsub(" ", "", strsplit(
      strsplit(iqtree_check_2, ":")[[1]][2], "\\("
    )[[1]][1])),
    "number_invariant_sites" = as.numeric(gsub(" ", "", strsplit(
      strsplit(iqtree_check_3, ":")[[1]][2], "\\("
    )[[1]][1])),
    "number_parsimony_informative_sites" = as.numeric(gsub(" ", "", strsplit(iqtree_check_4, ":")[[1]][2])),
    "number_distinct_site_patterns" = as.numeric(gsub(" ", "", strsplit(iqtree_check_5, ":")[[1]][2]))
  )
  iqtree_values <- c(
    iqtree_values,
    "percent_invariant_sites" = round(
      (iqtree_values[["number_invariant_sites"]] / iqtree_values[["num_sites"]] * 100),
      digits = 2),
    "percent_parsimony_informative_sites" = round(
      (iqtree_values[["number_parsimony_informative_sites"]] / iqtree_values[["num_sites"]] * 100),
      digits = 2)
  )
  return(iqtree_values)
}


check.tree.stats <- function(tree){
  # Calculate statistics on phylogenetic trees
  # Explanation of functions (from treestats doco):
  #   mntd: Per tip, evaluates the shortest distance to another tip, then takes the average across all tips
  #   mean_pair_dist: mean pairwise distance, using the fast algorithm by Constantinos, Sandel & Cheliotis (2012)
  #   var_pair_dist: After calculating all pairwise distances between all tips, this function takes the variance across these values.
  #   mean_branch_length: Mean branch length of a tree, including extinct branches
  #   var_branch_length: Variance of branch lengths of a tree, including extinct branches
  #   mean_branch_length_ext: Mean length of external branch lengths of a tree, e.g. of branches leading to a tip
  #   var_branch_length_ext: Variance of external branch lengths of a tree, e.g. of branches leading to a tip
  #   mean_branch_length_int: Mean length of internal branch lengths of a tree, e.g. of branches not leading to a tip
  #   var_branch_length_int: Variance of internal branch lengths of a tree, e.g. of branches not leading to a tip
  #   avg_vert_depth: The average vertex depth metric, measures the average path (in edges), between the tips and the root.
  #   max_width: Calculates the maximum width, this is calculated by first collecting the depth of each node and tip
  #               across the entire tree, where the depth represents the distance (in nodes) to the root. Then, the width
  #               represents the number of occurrences of each possible depth. The maximal width then returns the
  #               maximum number of such occurences
  #   max_depth: The maximum depth metric, measures the maximal path (in edges), between the tips and the root.
  #   max_ladder: Maximum ladder index (Higher values indicate more unbalanced tree)
  #   branching_times: Branching times of a tree
  #   crown_age: Crown age of a tree
  #   tree_height: Height of a tree (crown age)
  #   treeness: Calculates the fraction of tree length on internal branches, also known as treeness or stemmines
  tree_stats <- c(
    "mean_nearest_taxon_distance" = treestats::mntd(tree),
    "mean_pairwise_distance" = treestats::mean_pair_dist(tree),
    "mean_branch_length" = treestats::mean_branch_length(tree),
    "var_branch_length" = treestats::var_branch_length(tree),
    "mean_branch_length_external" = treestats::mean_branch_length_ext(tree),
    "var_branch_length_external" = treestats::var_branch_length_ext(tree),
    "mean_branch_length_internal" = treestats::mean_branch_length_int(tree),
    "var_branch_length_internal" = treestats::var_branch_length_int(tree),
    "average_vertex_depth" = treestats::avg_vert_depth(tree),
    "max_width" = treestats::max_width(tree),
    "max_depth" = treestats::max_depth(tree),
    "max_ladder" = treestats::max_ladder(tree),
    "max_branching_time" = max(treestats::branching_times(tree)),
    "crown_age" = treestats::crown_age(tree),
    "tree_height" = treestats::tree_height(tree),
    "treeness" = treestats::treeness(tree)
  )
  return(tree_stats)
}
