# /caitlinch/treelikeness-metrics/03_data_analysis.R
# Caitlin Cherryh 2022

# This program takes results from applying various treelikeness tests and performs data analysis



#### 1. Set parameters ####
# data_directory          <- Directory where alignments will be saved/treelikeness metrics will be run.
# output_directory        <- Directory for output of data analysis
# repo_directory          <- Location of caitlinch/treelikeness-metrics github repository (for access to functions).

data_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/01_results/"
output_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/02_data_analysis/"
repo_directory <- "/Users/caitlincherryh/Documents/Repositories/treelikeness_metrics/"


#### 2. Prepare analyses ####
# Open packages


# List all data files
data_files <- paste0(data_directory, list.files(data_directory))


# Create folder for plots
plot_directory <- paste0(output_directory, "plots/")
if (dir.exists(plot_directory) == FALSE){dir.create(plot_directory)}

#### 3. Analyse output from Experiment 1 ####
# Open data file from Experiment 1 as a dataframe
exp1_data_file <- grep("exp1", grep("treelikeness_metrics_results", data_files, value = TRUE), value = TRUE)
exp1_df <- read.csv(file = exp1_data_file, stringsAsFactors = FALSE)
