# /caitlinch/treelikeness-metrics/code/05_plotting.R
# Caitlin Cherryh 2022

# This program makes various additional or supplementary plots (main results are plotted in 04_data_analysis.R)



#### 1. Set parameters ####
# data_directory          <- Directory where alignments will be saved/treelikeness metrics will be run.
# output_directory        <- Directory for output of data analysis
# repo_directory          <- Location of caitlinch/treelikeness-metrics github repository (for access to functions)

data_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/01_results/"
output_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/02_data_analysis/"
repo_directory <- "/Users/caitlincherryh/Documents/Repositories/treelikeness-metrics/"



#### 2. Prepare analyses ####
# Source function files
source(paste0(repo_directory, "code/func_data_analysis.R"))

# Open packages
library(reshape2)
library(ggplot2)
library(patchwork)

# List all data files
data_files <- paste0(data_directory, list.files(data_directory))

# Create folder for plots
plot_directory <- paste0(output_directory, "plots/")
if (dir.exists(plot_directory) == FALSE){dir.create(plot_directory)}



#### 3. Plot timings for experiment 2 ####
# Identify the experiment 2 timings file
time_file <- grep("exp2", grep("test_statistic_timing", data_files, value = TRUE), value = TRUE)
# Open the timings file
time_df <- read.csv(time_file)
# Reshape into long format
long_df <- melt(time_df, id.vars = c("uid", "time_units"))
# Convert time to minutes
long_df$value <- round(long_df$value/60, digits = 2)
long_df$time_units <- "mins"
# Add fancy labels for facets
long_df$var_label <- factor(long_df$variable, 
                                 levels = c("tree_proportion", "Cunningham_test",
                                            "delta_plot", "likelihood_mapping",
                                            "ntlt", "q_residuals", 
                                            "scfs", "fast_tiger", 
                                            "total_time"), 
                                 ordered = TRUE, 
                                 labels = c(expression(atop("Tree","proportion")), expression(atop("Cunningham","metric")), 
                                            expression(paste(delta["q"])), expression(atop("Likelihood", "mapping")),
                                            expression("NTLT"), expression("Q-Residual"), 
                                            expression("sCF"), expression("Fast TIGER"),
                                            expression("Total")) )
# Create a discrete box plot where the tests are the x axis and time in minutes is on the y axis




