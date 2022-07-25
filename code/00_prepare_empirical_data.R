# /caitlinch/treelikeness-metrics/00_prepare_empirical_data.R
# Caitlin Cherryh 2022

# This program takes in a partition file and separates an alignment based on the charsets in the partition



#### 1. Set parameters ####
# alignment_file          <- nexus alignment file
# partition_file          <- nexus partition file 
# output_directory        <- Directory for output alignments
# repo_directory          <- Location of caitlinch/treelikeness-metrics github repository (for access to functions)

alignment_file <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/00_data_oaks2011/oaks2011_crocodylia_alignment.nex"
partition_file <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/00_data_oaks2011/oaks2011_crocodylia_mtDNA.nex"
output_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/00_data_oaks2011/mtDNA_gene_alignments/"
repo_directory <- "/Users/caitlincherryh/Documents/Repositories/treelikeness-metrics/"



#### 2. Prepare analyses ####
# Open libraries
