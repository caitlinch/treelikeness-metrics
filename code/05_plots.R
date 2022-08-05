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
# Add a column with the number of gene trees
time_df$num_genes <- 250
# Reshape into long format
long_df <- melt(time_df, id.vars = c("uid", "time_units", "num_genes"))
# Convert time to minutes
long_df$value <- round((long_df$value*(1/60)*(1/60)), digits = 2)
long_df$time_units <- "hours"
# Add fancy labels for facets
long_df$var_label <- factor(long_df$variable, 
                                 levels = c("likelihood_mapping", "scfs", "ntlt", "delta_plot", "q_residuals", "fast_tiger", "Cunningham_test", "tree_proportion", "total_time"), 
                                 ordered = TRUE, 
                                 labels = c("Likelihood\nmapping", "sCF", "Network\nTreelikeness\nTest", "Delta\nplot", "Q-Residual", "Fast\nTIGER", "Cunningham\nmetric", "Tree\nproportion", "Total\nruntime") )
# Create a discrete box plot where the tests are the x axis and time in minutes is on the y axis
p <- ggplot(data = long_df, aes(x = var_label, y = value)) +
  geom_boxplot(outlier.alpha = 0.6) +
  scale_x_discrete(name = "Test statistic") +
  scale_y_continuous(name = "Time to run test statistic (hours)", breaks = seq(0, 28, 4), labels = seq(0, 28, 4), minor_breaks = seq(0,28,2), limits = c(0,26)) +
  labs(title = "Experiment 2: Simulations with introgression (250 gene trees)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 17),
        axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12))
# Save plot
p_name <- paste0(plot_directory, "exp2_timing_boxplots_313obs_250genes.pdf")
ggsave(filename = p_name, plot = p, device = "pdf", height = 8, width = 8, units = "in")



