# caitlinch/treelikeness-metrics/code/05_figures.R
# Caitlin Cherryh 2023

# This program takes results from applying various treelikeness tests and performs data analysis/plotting


#### 1. Set parameters ####
# plot_directory           <- Directory for output of data analysis
# repo_directory          <- Location of caitlinch/treelikeness-metrics github repository (for access to functions)

plot_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/06_plots/"
repo_directory <- "/Users/caitlincherryh/Documents/Repositories/treelikeness-metrics/"



#### 2. Prepare analyses ####
# Open packages
library(ape)
library(phytools)
library(ggtree)
library(ggplot2)
library(patchwork)



#### 3. Plot figure showing ancient and recent recombination events ####
## Prepare tree for plotting
# Open a randomly generated tree. To make tree: rtree(10)
tree_text <- "(t6:0.07360411785,((t4:0.4142429121,(t5:0.9130682887,((t1:0.5073738024,(t10:0.8766066427,t8:0.2323632115):0.784319913):0.05937049002,(t2:0.2084415455,t9:0.7882211793):0.8298895105):0.3861808416):0.97176072):0.8059928031,(t3:0.04930912959,t7:0.06525251642):0.9080945631):0.3121802141);"
tree <- read.tree(text = tree_text)
# Root tree at taxa "t6"
tree <- root(tree, "t6", resolve.root = TRUE)
# Force tree to ultrametric
tree <- force.ultrametric(tree, method = "extend")
# Ladderize tree
tree <- ladderize(tree)
# Relabel tips from 1 - 10
tree$tip.label <- paste0("t", c("1", "4", "5", "8", "9", "10", "6", "7", "2", "3"))

## Plot a: Recent event
# Note: receptor = max(taxa_pair); donor = min(taxa_pair)
# For this example: "t5" -> "t7"
re_tree <- ggtree(tree, size = 1.5) +
  geom_rootedge(rootedge = 0.2, linewidth = 1.5) +
  geom_tiplab(offset = 0.1, size = 10) +
  xlim(-0.2,4.5) +
  geom_segment(aes(x = 3.77, y = 5.05, xend = 3.77, yend = 6.9), linewidth = 1.5, arrow = arrow(length = unit(0.7, "cm"), type = "closed"), color = "red")

## Plot b: Ancient event
# Note: receptor = max(taxa_pair); donor = min(taxa_pair)
# For this example: "t1" -> "t4"
ae_tree <- ggtree(tree, size = 1.5) +
  geom_rootedge(rootedge = 0.2, linewidth = 1.5) +
  geom_tiplab(offset = 0.1, size = 10) +
  xlim(-0.2,4.5) +
  geom_segment(aes(x = 0.73, y = 1.05, xend = 0.73, yend = 5.06), linewidth = 1.5, arrow = arrow(length = unit(0.7, "cm"), type = "closed"), color = "red")

## Combine tree plots with patchwork
quilt = re_tree + ae_tree + plot_annotation(tag_levels = 'a', tag_suffix = ")") & 
  theme(plot.tag = element_text(size = 40))
quilt_title <- paste0(plot_directory, "methods_introgression_events.")
ggsave(filename = paste0(quilt_title, "pdf"), plot = quilt, device = "pdf", width = 15, height = 6, units = "in")
ggsave(filename = paste0(quilt_title, "png"), plot = quilt, device = "png", width = 15, height = 6, units = "in")


