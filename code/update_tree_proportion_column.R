# Update the tree proportion column in experiment 1 and experiment 2

dir <- "~/Repositories/treelikeness-metrics/"

## For exp1:
# Get files
results_file <- paste0(dir, "output/exp1_treelikeness_metrics_collated_results.csv")
tp1_file <- paste0(dir, "output/exp1_TreeProportion_collated_results.csv")
tp2_file <- paste0(dir, "output/exp1_treeProportion_completeRuns_collated_results.csv")
# Open csvs
results_df <- read.csv(results_file)
tp1_df <- read.csv(tp1_file)
tp2_df <- read.csv(tp2_file)
# Check tp1 and tp2 have same tree proportion values
identical(tp1_df$uid,tp2_df$uid) # TRUE
identical(tp1_df$tree_proportion,tp2_df$tree_proportion) #TRUE
# tp1 and tp2 are identical. Ignore tp2 from now on.
# Compare results_df and tp1_df order of uids
identical(results_df$uid,tp1_df$uid) # FALSE
# Remove uids from results_df that are not needed
keep_df <- tp1_df[tp1_df$uid %in% results_df$uid, ]
# Check order of uids
identical(results_df$uid,tp1_df$uid) #FALSE
which((keep_df$uid == results_df$uid) == FALSE) # integer(0)
# Order of uids is the same
uid_check1 <- keep_df$uid
keep_df <- keep_df[match(results_df$uid, keep_df$uid),]
uid_check2 <- keep_df$uid
which((uid_check1 == uid_check2) == FALSE) # integer(0) (therefore no change in order of uids)
which(as.numeric(results_df$tree_proportion) != as.numeric(keep_df$tree_proportion)) # integer(0) (therefore all tp values identical for new function)
# Replace tree proportion values in the results csv
results_df$tree_proportion <- keep_df$tree_proportion
# Write out the results file
write.csv(results_df, file = results_file, row.names = FALSE)







