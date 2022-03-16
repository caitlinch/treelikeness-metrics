## Open packages


## Parameters for applying metrics
# local_directory         <- Directory where alignments will be saved/treelikeness metrics will be run.
# repo_directory          <- Location of caitlinch/treelikeness_metrics github repository (for access to functions).
# iqtree2_path            <- Path to IQ-Tree2.2-beta executable (this is the IQ-Tree2 release containing Alisim). 
# fast_TIGER_path         <- Path to fast TIGER executable.
# phylogemetric_path      <- Path to phylogemetric executable
# splitstree_path         <- Path to SplitsTree 4 version 4.17.2 or above
# netmake_path            <- Path to netmake executable within SPECTRE
# netme_path              <- Path to netme executable within SPECTRE

local_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/"
repo_directory <- "/Users/caitlincherryh/Documents/Repositories/treelikeness_metrics/"
iqtree2_path <- "iqtree2.2-beta"
fast_TIGER_path <- "/Users/caitlincherryh/Documents/Executables/fast_TIGER-0.0.2/DAAD_project/fast_TIGER"
phylogemetric_path <- "/Users/caitlincherryh/Documents/Executables/phylogemetric/phylogemetric_executable"
splitstree_path <- "/Applications/SplitsTree/SplitsTree.app/Contents/MacOS/JavaApplicationStub"
netmake_path <- "/Applications/Spectre.app/Contents/MacOS/netmake"
netme_path <- "/Applications/Spectre.app/Contents/MacOS/netme"


## Source functions from caitlinch/treelikeness_metrics
source(paste0(repo_directory, "func_metrics.R"))


## Find the three folders of simulated alignments
exp_folders <- paste0(local_directory, c("exp_1/", "exp_2/", "exp_3/"))

# For each experiment, get the list of folders within that experiment folder
exp1_runs <- paste0(exp_folders[1], list.files(exp_folders[1]), "/")
replicate_folder <- exp1_runs[1]
replicate_folder <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/exp_1/exp1_00010_0010_001_1/"

# test params
num_iqtree2_threads = "AUTO"
num_iqtree2_scf_quartets = 100
iqtree_substitution_model = "JC"
delta_plot_substitution_method = "JC69"
num_phylogemetric_threads = NA
sequence_format = "DNA"

treelikeness.metrics.simulations <- function(replicate_folder, iqtree2_path, splitstree_path, phylogemetric_path, fast_TIGER_path,
                                             num_iqtree2_threads = "AUTO", num_iqtree2_scf_quartets = 100, iqtree_substitution_model = "JC", 
                                             delta_plot_substitution_method = "JC69", num_phylogemetric_threads = NA, sequence_format = "DNA"){
  ## Function to take one alignment, apply all treelikeness metrics and return results in a dataframe
  
  # Get alignment file
  folder_files <- list.files(replicate_folder)
  alignment_path <- paste0(replicate_folder, grep("output_alignment\\.fa\\.", grep("output_alignment\\.fa", folder_files, value = TRUE), invert = TRUE, value = TRUE))
  df_name <- gsub("output_alignment.fa", "treelikeness_tests.csv", alignment_path)
  
  # Check whether dataframe .csv file already exists. If it does, import the dataframe. If it doesn't, make it.
  if (file.exists(df_name) == TRUE){
    ## Read in the results csv file
    results_df <- read.csv(df_name)
  } else if (file.exists(df_name) == FALSE){
    ## Apply all treelikeness test statistics to generate the results csv file
    
    # Determine the number of taxa (needed for number of quartets in likelihood mapping and sCFs)
    if (grepl("exp1", replicate_folder)){
      random_trees_file <- paste0(replicate_folder, grep("random_trees", folder_files, value = TRUE))
      random_trees <- read.tree(random_trees_file)
      n_tree_tips <- unique(Ntip(random_trees))[[1]]
    } else if ((grepl("exp2", replicate_folder))|(grepl("exp3", replicate_folder))){
      start_tree_file <- paste0(replicate_folder, grep("starting_tree", folder_files, value = TRUE))
      start_tree <- read.tree(start_tree_file)
      n_tree_tips <- Ntip(start_tree)
    }
    
    # Apply Likelihood mapping (Strimmer and von Haeseler 1997)
    lm <- likelihood.mapping(alignment_path, iqtree2_path, iqtree2_number_threads = num_iqtree2_threads, substitution_model = iqtree_substitution_model, 
                             number_of_taxa = n_tree_tips)
    
    # Apply Site concordance factors (Minh et. al. 2020)
    scfs <- scf(alignment_path, iqtree2_path, iqtree2_number_threads = num_iqtree2_threads, number_scf_quartets = num_iqtree2_scf_quartets, 
                substitution_model = iqtree_substitution_model, add.likelihood.map = FALSE, number_of_taxa = n_tree_tips)
    
    if (n_tree_tips <= 50){
      # Apply Network Treelikeness Test (Huson and Bryant 2006)
      ntlt <- network.treelikeness.test(alignment_path, splitstree_path, sequence_format)
    } else {
      # Do not apply Network Treelikeness Test - too many taxa
      ntlt <- NA
    }
    
    # Apply Delta plots (Holland et. al. 2002)
    mean_delta_plot_value <- mean.delta.plot.value(alignment_path, sequence_format, substitution_model = delta_plot_substitution_method)
    
    # Apply Q-residuals (Gray et. al. 2010)
    mean_q_residual <- q_residuals(alignment_path, phylogemetric_path, sequence_format, phylogemetric_number_of_threads = num_phylogemetric_threads)
    
    # Apply TIGER (Cummins and McInerney 2011)
    mean_tiger_value <- TIGER(alignment_path, fast_TIGER_path, sequence_format)
    
    # Apply Treeness test (Cavalli-Sforza and Piazza 1975)
    
    # Apply Cunningham test (Cunningham 1975)
    
    # Apply tree proportion (new test)
    
    ## Assemble results into a dataframe and save
    results_vec <- c(lm, scfs$mean_scf, scfs$median_scf, min(scfs$all_scfs), max(scfs$all_scfs), ntlt, mean_delta_plot_value, mean_q_residual, mean_tiger_value)
    results_df <- as.data.frame(matrix(data = results_vec, nrow = 1, ncol = length(results_vec), byrow = TRUE))
    names_vec <- c("LM_num_resolved_quartets", "LM_num_partly_resolved_quartets", "LM_num_unresolved_quartets", "LM_total_num_quartets", "LM_proportion_resolved_quartets",
                   "sCF_mean", "sCF_median", "sCF_min", "sCF_max", "NetworkTreelikenessTest", "NTLT_AllSplitsCompatible", "NTLT_proportion_compatible_split_comparisons",
                   "mean_delta_plot_value", "mean_Q_residual", "mean_TIGER_value")
    names(results_df) <- names_vec
    write.csv(results_df, file = df_name, row.names = FALSE)
  }
  
  # Return results_df so it can be collated using lapply
  return(results_df)
}

