# /caitlinch/treelikeness-metrics/code/00_prepare_empirical_data.R
# Caitlin Cherryh 2022

# This program takes in a partition file and separates an alignment based on the charsets in the partition

# Paper:
#   Oaks, J.R. (2011), A time-calibrated species tree of Crocodylia reveals a recent radiation of the true crocodiles. Evolution, 65: 3285-3297. https://doi.org/10.1111/j.1558-5646.2011.01373.x

# Data available from: 
#   Oaks, Jamie R (2011), Data from: A time-calibrated species tree of Crocodylia reveals a recent radiation of the true crocodiles, Dryad, Dataset, https://doi.org/10.5061/dryad.5k9s0

#### 1. Set parameters ####
# alignment_file          <- nexus alignment file, downloaded from the DataDryad above
# output_directory        <- Directory for saving output alignments
# repo_directory          <- Location of caitlinch/treelikeness-metrics github repository

alignment_file <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/00_data_oaks2011/oaks2011_crocodylia_alignment.nex"
output_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/00_data_oaks2011/"
repo_directory <- "/Users/caitlincherryh/Documents/Repositories/treelikeness-metrics/"



#### 2. Prepare analyses ####
# Open libraries
library(phylotools)
library(phangorn)

# Convert nexus to fasta
fasta_file <- gsub("\\.nex", "\\.fa", alignment_file)
# Open alignment file as matrix
dna_nex <- read.nexus.data(alignment_file)
write.dna(dna_nex, file = fasta_file, format = "fasta", colsep = "")
# Open fasta file
dna_mat <- as.matrix(read.dna(file = fasta_file, format = "fasta"))

# Find partition files
all_empirical_files <- list.files(paste0(repo_directory, "empirical_analysis_files/"))
oaks_files <- paste0(repo_directory, "empirical_analysis_files/", grep("oaks2011", all_empirical_files, value = TRUE))
partition_files <- grep("partition", oaks_files, value = TRUE)


#### 3. Save each gene individually ####
charset_type = c("mtDNA", "nDNA")
for (c in charset_type){
  # Make new output directory
  gene_directory <- paste0(output_directory, c, "_gene_alignments/")
  if (dir.exists(gene_directory) == FALSE){dir.create(gene_directory)}
  
  # Select partition file
  c_partition_file <- grep(c, partition_files, value = TRUE)
  # Open partition file as lines
  partitions <- readLines(c_partition_file)
  # Extract only lines that contain a charset and remove extra formatting
  partitions <- grep("charset", partitions, value = TRUE)
  partitions <- gsub("\tcharset ", "", partitions)
  partitions <- gsub(";", "", partitions)
  
  ## Save each gene
  print('Saving individual genes:')
  all_gene_lengths <- c()
  all_gene_names <- c()
  for (i in 1:length(partitions)){
    # Get the ith gene from the list of genes
    gene_line <- partitions[i]
    
    # Split the line to determine the gene name
    gene_line <- gsub(";", "", gene_line)
    gene_line_split <- strsplit(gene_line, "=")[[1]]
    gene_name <- gsub(" ", "", gene_line_split[1])
    print(gene_name)
    
    # Check whether this line is for a single codon position or a whole gene
    get_codon_position <- grepl("\\\\3", gene_line)
    
    # Split the line to determine the gene start and end position
    gene_range <- gsub(" ", "", gene_line_split[2])
    gene_range_split <- strsplit(gene_range, "-")[[1]]
    if (length(gene_range_split) == 2){
      gene_start <- as.numeric(gene_range_split[[1]])
      if (get_codon_position == FALSE){
        # If getting whole gene, just take end of range
        gene_end <- as.numeric(gene_range_split[[2]])
        # Get the whole range of numbers as the gene_positions to extract
        gene_positions <- gene_start:gene_end
      } else if (get_codon_position == TRUE){
        # If getting only one codon position, get only last number (do not include the "\\3")
        gene_end <- as.numeric(gsub("\\\\3", "", gene_range_split[[2]]))
        # Get every third number starting from the gene_start and moving towards the gene_end
        gene_positions <- seq(gene_start, gene_end, 3)
      }
    } else if (length(gene_range_split) > 2){
      # If there are multiple components to the partition, extract all the sites for inclusion
      all_ranges <- strsplit(gene_line_split[2], " ")[[1]]
      all_ranges <- all_ranges[which(all_ranges != "")]
      gene_positions <- c()
      for (ind in 1:length(all_ranges)){
        j <- all_ranges[ind]
        j_get_codon_position <- grepl("\\\\3", j)
        if (j_get_codon_position == FALSE){
          j_split <- strsplit(j, "-")[[1]]
          j_start <- as.numeric(j_split[[1]])
          j_end <- as.numeric(j_split[[2]])
          j_positions <- c(j_start:j_end)
          gene_positions <- c(gene_positions, j_positions)
        } else if  (j_get_codon_position == TRUE){
          j_2 <- gsub("\\\\3", "", j)
          j_split <- strsplit(j_2, "-")[[1]]
          j_start <- as.numeric(j_split[[1]])
          j_end <- as.numeric(j_split[[2]])
          j_positions <- seq(j_start, j_end, 3)
          gene_positions <- c(gene_positions, j_positions)
        }
      }
    }
    
    # Length of gene will just be the number of sites included in the gene_positions
    gene_length <- length(gene_positions)
    
    # Save name and length of gene
    all_gene_lengths <- c(all_gene_lengths, gene_length)
    all_gene_names <- c(all_gene_names, gene_name)
    
    # Sort gene positions into numerical order
    # Subset the supermatrix to get the sites (columns) for this gene
    sorted_gene_positions <- sort(gene_positions, decreasing = FALSE)
    gene_mat <- dna_mat[, c(sorted_gene_positions)]
    
    # Assemble file name for gene
    gene_file <- paste0(gene_directory, gene_name, ".fa")
    # Write gene to file
    write.FASTA(gene_mat, file = gene_file)
  }
  
  # Create dataframe with gene name and lengths
  df <- data.frame(dataset = "Oaks2011", gene_name = all_gene_names, gene_length = all_gene_lengths)
  write.csv(df, file = paste0(output_directory, "Oaks2011_", c,"_gene_lengths.csv"), row.names = FALSE)
}



#### 3. Create subsets of taxa for each gene  ####
# Need to separate each alignment into the three subsets for analysis

# This subset contains every individual from the alignment
all_individuals <- c("LSUMZ_H18733", 
                     "LSUMZ_H7868", "LSUMZ_H21699", "LSUMZ_H21700",
                     "LSUMZ_H13961", "LSUMZ_H13962", "LSUMZ_H13964", "LSUMZ_H21701", "LSUMZ_H21702", 
                     "LSUMZ_H21705", "LSUMZ_H21706",
                     "LSUMZ_H21707", 
                     "LSUMZ_H21751", "LSUMZ_H21752", "LSUMZ_H21753",
                     "LSUMZ_H6998", "LSUMZ_H21761", 
                     "LSUMZ_H6420", "LSUMZ_H7873",
                     "LSUMZ_H21748", 
                     "LSUMZ_H21763", "LSUMZ_H21764", "LSUMZ_H21765", 
                     "LSUMZ_H6760", "LSUMZ_H6982", "LSUMZ_H21708", "LSUMZ_H21709", "LSUMZ_H21710", "LSUMZ_H21711", "LSUMZ_H21712", "LSUMZ_H21713", "LSUMZ_H21714", "LSUMZ_H21715",
                     "LSUMZ_H6976", "LSUMZ_H21718", "LSUMZ_H21719", "LSUMZ_H21720", 
                     "LSUMZ_H20683", "LSUMZ_H20684", "LSUMZ_H20685", "LSUMZ_H20686", "LSUMZ_H21724", 
                     "LSUMZ_H7070", "LSUMZ_H21725", "LSUMZ_H21726",
                     "LSUMZ_H6903", "LSUMZ_H21727", "LSUMZ_H21729", "LSUMZ_H21730",
                     "LSUMZ_H21731", "LSUMZ_H21733", "LSUMZ_H21734", "LSUMZ_H21735", "LSUMZ_H21736", "LSUMZ_H21737", "LSUMZ_H21738", "LSUMZ_H21739",
                     "LSUMZ_H21766", "LSUMZ_H21768", "LSUMZ_H21769", "LSUMZ_H21831", "LSUMZ_H21771", "LSUMZ_H21815", "LSUMZ_H21872",
                     "LSUMZ_H6995", "LSUMZ_H7071", 
                     "LSUMZ_H21741", "LSUMZ_H21742",
                     "LSUMZ_H6758", "LSUMZ_H6984",
                     "LSUMZ_H21745", "LSUMZ_H21746", "LSUMZ_H21747", 
                     "LSUMZ_H6978", "LSUMZ_H6985", 
                     "LSUMZ_H21755", "LSUMZ_H21756", "LSUMZ_H6990", "LSUMZ_H6992")
all_species <- c("Alligator mississippiensis", 
                 rep("Alligator sinensis", 3),
                 rep("Caiman crocodilus", 5),
                 rep("Caiman latirostris", 2),
                 "Caiman yacare", 
                 rep("Melanosuchus niger", 3),
                 rep("Paleosuchus palpebrosus", 2),
                 rep("Paleosuchus trigonatus", 2), 
                 "Gavialis gangeticus", 
                 rep("Tomistoma schlegelii", 3),
                 rep("Crocodylus acutus", 10),
                 rep("Mecistops cataphractus", 4), 
                 rep("Crocodylus intermedius", 5),
                 rep("Crocodylus johnstoni", 3),
                 rep("Crocodylus moreletii", 4),
                 rep("Crocodylus niloticus", 8),
                 rep("Crocodylus mindorensis", 7),
                 rep("Crocodylus novaeguineae", 2),
                 rep("Crocodylus palustris", 2),
                 rep("Crocodylus porosus", 2),
                 rep("Crocodylus rhombifer", 3),
                 rep("Crocodylus siamensis", 2),
                 rep("Osteolaemus tetraspis", 4) )
# This subset contains one individual per species (the individual present within the alignment with the lowest LSUMZ ID number)
one_individuals <- c("LSUMZ_H18733", "LSUMZ_H7868", "LSUMZ_H13961", "LSUMZ_H21705",
                     "LSUMZ_H21707", "LSUMZ_H21751", "LSUMZ_H6998", "LSUMZ_H6420",
                     "LSUMZ_H21748", "LSUMZ_H21763", "LSUMZ_H6760", "LSUMZ_H6976",
                     "LSUMZ_H20683", "LSUMZ_H7070", "LSUMZ_H6903", "LSUMZ_H21731",
                     "LSUMZ_H21766", "LSUMZ_H6995", "LSUMZ_H21741", "LSUMZ_H6758", 
                     "LSUMZ_H21745", "LSUMZ_H6978", "LSUMZ_H6990")
one_species <- c("Alligator mississippiensis", "Alligator sinensis", "Caiman crocodilus", "Caiman latirostris",
                 "Caiman yacare", "Melanosuchus niger", "Paleosuchus palpebrosus", "Paleosuchus trigonatus",
                 "Gavialis gangeticus", "Tomistoma schlegelii", "Crocodylus acutus", "Mecistops cataphractus",
                 "Crocodylus intermedius", "Crocodylus johnstoni", "Crocodylus moreletii", "Crocodylus niloticus",
                 "Crocodylus mindorensis", "Crocodylus novaeguineae", "Crocodylus palustris", "Crocodylus porosus",
                 "Crocodylus rhombifer", "Crocodylus siamensis", "Osteolaemus tetraspis")
# This subset contains one individual per species (the individual present within the alignment with the lowest LSUMZ ID number) for a single clade within the main tree
clade_individuals <- c("LSUMZ_H18733", "LSUMZ_H7868", "LSUMZ_H13961", "LSUMZ_H21705",
                       "LSUMZ_H21707", "LSUMZ_H21751", "LSUMZ_H6998", "LSUMZ_H6420")
clade_species <- c("Alligator mississippiensis", "Alligator sinensis", "Caiman crocodilus", "Caiman latirostris",
                   "Caiman yacare", "Melanosuchus niger", "Paleosuchus palpebrosus", "Paleosuchus trigonatus")

# Create a list of all the information for each subset
subset_list = list(all = list(individuals = all_individuals, species = all_species, num_taxa = length(all_individuals), id = paste0(length(all_individuals), "taxa")),
                   species = list(individuals = one_individuals, species = one_species, num_taxa = length(one_individuals), id = paste0(length(one_individuals), "taxa")),
                   clade = list(individuals = clade_individuals, species = clade_species, num_taxa = length(clade_individuals), id = paste0(length(clade_individuals), "taxa")) )
keys = c("all", "species", "clade") 

# Iterate through each of the three subsets
for (k in keys){
  # Iterate through each of the charsets
  charset_type = c("mtDNA", "nDNA")
  for (c in charset_type){
    # Make new output directory
    gene_directory <- paste0(output_directory, c, "_gene_alignments/")
    if (dir.exists(gene_directory) == FALSE){dir.create(gene_directory)}
    
    # List all files in that directory
    all_genes <- paste0(gene_directory, list.files(gene_directory))
    # Remove any genes with any of the ids
    all_genes <- grep(subset_list[[1]]$id, grep(subset_list[[2]]$id, grep(subset_list[[3]]$id, all_genes, invert = TRUE, value = TRUE), invert = TRUE, value = TRUE), invert = TRUE, value = TRUE)
    
    # Get the list for k
    k_list <- subset_list[[k]]
    # Get the id for k
    k_id <- k_list[["id"]]
    # Get the list of taxa for k
    k_individuals <- k_list[["individuals"]]
    # Take the whole list of taxa and return those that are NOT in k_individuals but ARE in the whole list
    remove_taxa <- setdiff(subset_list$all$individuals, k_individuals)
    
    # For each gene, open that gene and remove all unnecessary taxa
    lapply(all_genes, remove.unnecessary.taxa, taxa_to_remove = remove_taxa, id = k_id)
  }
}
