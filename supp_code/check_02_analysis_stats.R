## Check all of the alignments in the simulations (and the trees used to
#   simulate them) to determine how often we have this problem.
#   * To do this, you should re-run code from across your alignments, then examine
#     them programmatically (e.g. by plotting out tree lengths and actual depths,
#     by getting alignment statistics like pairwise identity etc).
#   * You should also examine (by eye) alignments from across all the simulations
#     to see if there are other issues we aren't aware of yet.

## 01. Prepare input parameters
num_cores <- 10
out_dir <- "/mnt/data/dayhoff/home/u5348329/treelikeness_supp/suppA/"

## 02. Prepare libraries and functions
library(parallel)
source("/mnt/data/dayhoff/home/u5348329/treelikeness_supp/supp_code/check_funcs.R")


## 03. Collect expA1 simulation statistics
expA1_dir <- "/mnt/data/dayhoff/home/u5348329/treelikeness_supp/suppA/exp_A1"
expA1_dirs <- paste0(list.dirs(expA1_dir), "/")
expA1_dirs <- expA1_dirs[2:length(expA1_dirs)]
mclapply(
  expA1_dirs,
  check.expA1.analysis.stats,
  mc.cores = num_cores
)
expA1_check_files <- paste0(expA1_dir,
                            "/",
                            grep(
                              ".check_02_analysis_stats.csv",
                              list.files(expA1_dir, recursive = TRUE),
                              value = TRUE
                            ))
expA1_collated <- as.data.frame(do.call(rbind, lapply(expA1_check_files, read.csv)))
write.csv(expA1_collated,
          file = paste0(out_dir, "expA1_check_02_analysis_stats_collated.csv"),
          row.names = FALSE)


## 04. Collect expA2 simulation statistics
expA2_dir <- "/mnt/data/dayhoff/home/u5348329/treelikeness_supp/suppA/exp_A2"
expA2_dirs <- paste0(list.dirs(expA2_dir), "/")
expA2_dirs <- expA2_dirs[2:length(expA2_dirs)]
mclapply(
  expA2_dirs,
  check.expA2.analysis.stats,
  mc.cores = num_cores
)
expA2_check_files <- paste0(expA2_dir,
                            "/",
                            grep(
                              ".check_02_analysis_stats.csv",
                              list.files(expA2_dir, recursive = TRUE),
                              value = TRUE
                            ))
expA2_collated <- as.data.frame(do.call(rbind, lapply(expA2_check_files, read.csv)))
write.csv(expA2_collated,
          file = paste0(out_dir, "expA2_check_02_analysis_stats_collated.csv"),
          row.names = FALSE)

