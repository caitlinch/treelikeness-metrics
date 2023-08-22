# caitlinch/treelikeness_metrics/code/func_parametric_bootstraps.R
# Caitlin Cherryh 2023

# This file contains functions for generating parametric bootstraps from an empirical alignment

library(ape)




#### Gap adding functions ####
copy.alignment.gaps <- function(alignment_path, template_alignment_path){
  ## Function to take a template alignment and copy the gaps onto another alignment
  
  # Open both the alignment and template alignment
  al <- as.character(read.FASTA(alignment_path, type = "AA"))
  template_al <- as.character(read.FASTA(template_alignment_path, type = "AA"))
  # Get the names of all the sequences
  seq_names <- names(al)
  # Iterate through each row and add gaps into the simulated alignment path
  for (seq_name in seq_names){
    original_seq <- template_al[[seq_name]] # get the original empirical sequence
    new_seq <- al[[seq_name]] # get the new simulated sequence that has the same name
    gap_inds <- which(original_seq == "-") # find out which sites are a gap in the original alignment
    unknown_inds <- which(original_seq == "?") # find out which sites are unknown in the original alignment
    new_seq[gap_inds] <- "-" # add the gaps into the simulated alignment
    new_seq[unknown_inds] <- "?" # add the unknowns into the simulated alignment
    al[[seq_name]] <- new_seq
  }
  # Convert alignment back to AAbin class
  al <- as.AAbin(al)
  # Save updated alignment to the alignment path
  write.FASTA(al, file = alignment_path, header = NULL, append = FALSE)
  # Return alignment path
  return(alignment_path)
}



#### Model extraction functions ####
identify.best.model <- function(iqtree_file){
  ## Collate best model and rate from iqtree file
  
  # Extract model
  model <- extract.best.model(iqtree_file)
  # Extract rates
  rates <- extract.rates(iqtree_file)
  # Check whether both exist
  if (is.na(model) == FALSE & is.na(rates) == FALSE){
    ## Both model and rates exist
    # Break up model and remove the rate category argument
    model_split <- strsplit(model, "\\+")[[1]]
    other_model_chunks <- grep("R1|R2|R3|R4|R5|R6|R7|R8|R9|R10", model_split, value = T, invert = T)
    rate_model_chunks <- grep("R1|R2|R3|R4|R5|R6|R7|R8|R9|R10", model_split, value = T, invert = F)
    # Collate all bits of the model into one string
    best_model <- paste0("'", paste(other_model_chunks, collapse = "+"), "+", rate_model_chunks, "{", rates, "}", "'")
  } else if (is.na(model) == FALSE & is.na(rates) == TRUE){
    ## Only model exists (no rates)
    best_model <- model
  } else {
    ## In any other case, return NA
    best_model <- NA
  }
  
  # Return the best model
  return(best_model)
}



extract.best.model <- function(iqtree_file){
  ### Extracts the best model of sequence evolution given a .iqtree file
  
  # Check if the iqtree file exists
  if (file.exists(iqtree_file) == TRUE){
    # If the iqtree_file does exist:
    ## Open the .iqtree file:
    iq_lines <- readLines(iqtree_file)
    
    ## Check for a ModelFinder section:
    # Determine whether there is a ModelFinder section
    mf_ind <- grep("ModelFinder", iq_lines)
    # Determine whether there is a line detailing the best model
    bm_ind <- grep("Best-fit model according to", iq_lines)
    
    ## Check for a Substitution Process section:
    # Determine the starting line of this section
    sp_ind <- grep("SUBSTITUTION PROCESS", iq_lines)
    # Determine the line detailing the model used
    mos_ind <- grep("Model of substitution", iq_lines)
    
    ## Extract the best fit model from the .iqtree file:
    if ((identical(mf_ind, integer(0)) == FALSE) & (identical(bm_ind, integer(0)) == FALSE)){
      # If ModelFinder was run, extract the best model from the ModelFinder section of the .iqtree file
      # Extract the line containing the best fit model
      m_line <- iq_lines[bm_ind]
    } else if ((identical(sp_ind, integer(0)) == FALSE) & (identical(mos_ind, integer(0)) == FALSE)) {
      # If there is no ModelFinder section, extract the model used from the substitution process section
      m_line <- iq_lines[mos_ind]
    } else {
      m_line <- "NA:NA"
    }
    
    ## Format the model nicely for output: 
    # Split the line at the colon into two parts
    m_line_split <- strsplit(m_line, ":")[[1]]
    # If the best model is a single model, the length of m_line_split will be 2
    #     One section for the explanatory text and one for the model
    # If the best model is a partition model, it will have more than two sections when split by colons
    # Extract the second part of the line onwards (contains the best fit model)
    best_model <- m_line_split[2:length(m_line_split)]
    # If best_model is longer than 1, paste it together again using colons
    if (length(best_model) >1){
      best_model <- paste(best_model, collapse = ":")
    }
    # Remove any white space from the best model
    best_model <- gsub(" ", "", best_model)
  } else if (file.exists(iqtree_file) == FALSE){
    # If the iqtree_file doesn't exist, return NA
    best_model = NA
  } # end if (file.exists(iqtree_file) == TRUE){
  
  # Return the best model from the iqtree_file (if the file exists)
  return(best_model)
}



extract.rates <- function(iqtree_file){
  # Function to extract the rate parameters of a model from in IQ-Tree
  
  # Check if the iqtree file exists
  if (file.exists(iqtree_file) == TRUE){
    # If the iqtree_file does exist:
    ## Open the .iqtree file:
    iq_lines <- readLines(iqtree_file)
    
    # Check for +R parameters 
    rate_ind <- grep("Site proportion and rates\\:", iq_lines)
    # Check if the rate_ind is present (if yes, that means there's a line containing the rates and weights)
    if (identical(rate_ind,integer(0)) == FALSE & class(rate_ind) == "integer"){
      ## Rates (+R)
      # The site proportion and weights values are present
      # Extract the rate weights and parameters from the iqtree file
      rates_line <- iq_lines[rate_ind]
      # Remove text from the beginning of the line
      rates_raw <- strsplit(rates_line, "\\:")[[1]][2]
      # Split up by the spaces
      split_rates <- strsplit(rates_raw, " ")[[1]]
      # Remove any entries that are empty characters (i.e. "")
      split_rates <- split_rates[split_rates != ""]
      # Split rates again by the commas (",")
      split2_rates <- unlist(strsplit(split_rates, ","))
      # Remove brackets from values
      rate_vals <- unlist(strsplit(unlist(strsplit(split2_rates, "\\(")), "\\)"))
      # Make output (for attaching into IQ-Tree2 command line)
      # Order is: proportion_1, rate_1, proportion_2, rate_2, ....., proportion_n, rate_n
      rates_op <- paste(rate_vals, collapse = ",")
    } else {
      # The site proportion and weights values are missing
      # Return NA for this file
      rates_op <- NA
    } # end if (identical(rate_ind,logical(0)) == FALSE & class(rate_ind) == "integer")
  } # end if (file.exists(iqtree_file) == TRUE)
  
  # Return the output
  return(rates_op)
}



get.simulation.parameters <- function(dotiqtree_file){
  ### Extracts all possible parameters from a .iqtree file
  
  # read in the IQ-TREE file to get substitution model and parameters
  iq_file <- readLines(dotiqtree_file)
  # extract the file name
  ind      <- grep("Input file name:",iq_file)
  op1      <- substr(iq_file[[ind]],18,nchar(iq_file[[ind]]))
  # extract the number of taxa and extract the length of the alignment
  ind         <- grep("Input data:",iq_file)
  input_str   <- iq_file[[ind]] # get the line that contains this info
  input_ls    <- strsplit(input_str," ")
  op2         <- input_ls[[1]][3] # extract number of sequences (number of taxa)
  op3         <- input_ls[[1]][6] # extract number of sites 
  # Extract the model of substitution (same for amino-acid and nucleotide files)
  ind         <- grep("Model of substitution: ",iq_file)
  sub_str     <- iq_file[[ind]] # get the line that contains this info
  sub_ls      <- strsplit(sub_str," ")
  op4         <- sub_ls[[1]][4]
  # Extract information about the sequence alignment
  ind <- grep("Number of constant sites:",iq_file)
  num_lines <- iq_file[c(ind:(ind+3))] # take the four lines listing the number of different kinds of sites
  num_split <- unlist(strsplit(num_lines,":")) # Split the lines at the colon
  num_names <- num_split[c(TRUE,FALSE)] # before the colon = name
  num_vals <- num_split[c(FALSE,TRUE)] # after the colon = value
  num_vals_regx <- regmatches(num_vals, gregexpr("[[:digit:]]+", num_vals)) # extract all the numbers after the colon
  # four strings = four lists of numbers: take first value from each list to get number of sites
  num_vals <- c(num_vals_regx[[1]][1],num_vals_regx[[2]][1],num_vals_regx[[3]][1],num_vals_regx[[4]][1]) 
  
  # check the type of sites - amino acid or DNA
  # if the input is DNA (nucleotide sites), gather that information
  if (input_ls[[1]][7]=="nucleotide"){
    # Extract the rate parameters
    rate1 <- as.numeric(strsplit(iq_file[[grep("A-C",iq_file)]],":")[[1]][2]) # A-C rate (same as code above, but combined 4 lines into 1 line)
    rate2 <- as.numeric(strsplit(iq_file[[grep("A-G",iq_file)]],":")[[1]][2]) # A-G rate
    rate3 <- as.numeric(strsplit(iq_file[[grep("A-T",iq_file)]],":")[[1]][2]) # A-T rate
    rate4 <- as.numeric(strsplit(iq_file[[grep("C-G",iq_file)]],":")[[1]][2]) # C-G rate
    rate5 <- as.numeric(strsplit(iq_file[[grep("C-T",iq_file)]],":")[[1]][2]) # C-T rate
    rate6 <- as.numeric(strsplit(iq_file[[grep("G-T",iq_file)]],":")[[1]][2]) # G-T rate
    
    # Extract the state frequencies
    state_freq_line <- iq_file[[grep("State frequencies",iq_file)]]
    if (state_freq_line == "State frequencies: (equal frequencies)"){
      # If the state frequencies are all equal, assign them all to 0.25 (1/4)
      sf1 <- 0.25 # pi(A) - A freq.
      sf2 <- 0.25 # pi(C) - C freq.
      sf3 <- 0.25 # pi(G) - G freq.
      sf4 <- 0.25 # pi(T) - T freq.
    } else {
      # If the state frequencies are not all equal, extract what they are
      sf1 <- as.numeric(strsplit(iq_file[[grep("pi\\(A\\)",iq_file)]],"=")[[1]][2]) # pi(A) - A freq. Remember to double backslash to escape before brackets
      sf2 <- as.numeric(strsplit(iq_file[[grep("pi\\(C\\)",iq_file)]],"=")[[1]][2]) # pi(C) - C freq.
      sf3 <- as.numeric(strsplit(iq_file[[grep("pi\\(G\\)",iq_file)]],"=")[[1]][2]) # pi(G) - G freq.
      sf4 <- as.numeric(strsplit(iq_file[[grep("pi\\(T\\)",iq_file)]],"=")[[1]][2]) # pi(T) - T freq.
    }
    
    # Extract model of rate heterogeneity
    mrh1      <- strsplit(iq_file[[grep("Model of rate heterogeneity:",iq_file)]],":")[[1]][2] # Extract model of rate heterogeneity 
    mrh2      <- as.numeric(strsplit(iq_file[[(grep("Model of rate heterogeneity:",iq_file)+1)]],":")[[1]][2]) # Line after the "model of rate heterogeneity" varies - extract it regardless of what it is 
    mrh2_name <- strsplit(iq_file[[(grep("Model of rate heterogeneity:",iq_file)+1)]],":")[[1]][1] # As the line varies, extract the name for the output dataframe
    mrh2_name <- gsub(" ","_",mrh2_name) # change the name to be easy to parse
    
    # make a list of the rows for the output dataframe
    names <- c("file_name","sequence_type","n_taxa","n_sites",num_names,"substitution_model", "A-C_rate","A-G_rate","A-T_rate","C-G_rate","C-T_rate","G-T_rate","A_freq","C_freq",
               "G_freq","T_freq","model_of_rate_heterogeneity","model_of_rate_heterogeneity_line2_name","model_of_rate_heterogeneity_line2_value")
    # Make a list of the output rows for the output dataframe
    op <- c(op1,"DNA",op2,op3,num_vals,op4,rate1,rate2,rate3,rate4,rate5,rate6,sf1,sf2,sf3,sf4,mrh1,mrh2_name,mrh2)
    # Create the output dataframe
    op_df <- data.frame(names,op, stringsAsFactors = FALSE)
    # Name the columns
    names(op_df) <- c("parameter","value")
    
    # Create the rate matrix Q
    Q_start <- grep("Rate matrix Q:",iq_file)+2
    Q_end   <- Q_start+3
    # Create the columns
    c1 <- c("A","C","G","T")
    c2 <- c()
    c3 <- c()
    c4 <- c()
    c5 <- c()
    # For each row in the iqtree file rate matrix
    for (i in Q_start:Q_end){
      # Split the row
      row <- strsplit(iq_file[[i]]," ")[[1]]
      row <- row[str_detect(row,"([0-9])")] # take only the numeric elements of the vector
      # Add the resulting values to the relevant columns
      c2 <- c(c2,as.numeric(row[1])) # convert to numeric so can use the numbers more easily later
      c3 <- c(c3,as.numeric(row[2]))
      c4 <- c(c4,as.numeric(row[3]))
      c5 <- c(c5,as.numeric(row[4]))
    }
    # Create a dataframe of the rate matrix Q
    q_df <- data.frame(c1,c2,c3,c4,c5, stringsAsFactors = FALSE)
    #Rename the columns
    names(q_df) <- c("nucleotide","A","C","G","T")
    
    # Check if the model for rate heterogeneity is uniform
    mrh1_check <- gsub(" ","",mrh1)
    if (mrh1_check=="Uniform"){
      # If the model for rate heterogeneity is uniform, don't need to create a matrix for discrete gamma rate categories
      g_df <- "Uniform"
    } else {
      # If the model isn't uniform, need to create a matrix to collect and store the gamme category information
      #Create the matrix for discrete gamma categories
      g_start <- grep(" Category",iq_file)+1 # get the index for the first line of the gamma categories matrix
      empty   <- which(iq_file=="") # get indexes of all empty lines
      empty   <- empty[empty>g_start] # get empty lines above gamma categories matrix
      g_end   <- empty[1]-1 # get end index for gamma categories matrix (one less than next empty line)
      end_line <- iq_file[g_end]
      # if the end isn't an empty line, subtract one from the end count 
      # to exclude lines like "Relative rates are computed as MEAN of the portion of the Gamma distribution falling in the category."
      # to see if this is what's happening, check whether the line starts with a numeric section (i.e. a category for the gamma rate)
      check_line <- length(strsplit(strsplit(end_line, "        " )[[1]][1]," ")[[1]])
      if (check_line > 3){
        # If the check_line is longer than 3 characters, it won't be a group for the gamma categories but an instruction
        # Instructions can be excluded from the gamma matrix (but categories can't)
        g_end = g_end - 1
      }
      # Start collecting info for the matrix
      g1 <- c() # initialise columns to store data in
      g2 <- c()
      g3 <- c()
      # Iterate through rows in gamma matrix
      for (i in g_start:g_end){
        row <- strsplit(iq_file[[i]]," ")[[1]] # split the rows on the (large amount of) " "'s in the middle
        row <- row[str_detect(row,"([0-9])")] # take only the numeric elements of the vector
        g1 <- c(g1,as.numeric(row[1])) # add the values to the columns
        g2 <- c(g2,as.numeric(row[2]))
        g3 <- c(g3,as.numeric(row[3]))
      }
      g_df <- data.frame(g1,g2,g3, stringsAsFactors = FALSE) # create a dataframe of the information
      names(g_df) <- c("category","relative_rate","proportion") # name the columns
    }
    
    # Create a list of the three dataframes
    # This will be the output 
    params <- list(op_df,g_df,q_df)
    # Name the parameters so they're easy to access once you've outputted the data
    names(params) <- c("parameters","gamma_categories","Q_rate_matrix")
    
  } else if (input_ls[[1]][7]=="amino-acid"){ # alternatively if the data is amino acid sites
    # Extract model of rate heterogeneity
    mrh1      <- strsplit(iq_file[[grep("Model of rate heterogeneity:",iq_file)]],":")[[1]][2] # Extract model of rate heterogeneity 
    mrh2      <- as.numeric(strsplit(iq_file[[(grep("Model of rate heterogeneity:",iq_file)+1)]],":")[[1]][2]) # Line after the "model of rate heterogeneity" varies - extract it regardless of what it is 
    mrh2_name <- strsplit(iq_file[[(grep("Model of rate heterogeneity:",iq_file)+1)]],":")[[1]][1] # As the line varies, extract the name for the output dataframe
    # Extract state frequencies
    sf1      <- strsplit(iq_file[[grep("State frequencies:",iq_file)]],":")[[1]][2]
    
    # make a list of the rows for the output dataframe
    names <- c("file_name","sequence_type","n_taxa","n_sites",num_names,"substitution_model","model_of_rate_heterogeneity","model_of_rate_heterogeneity_line2_name",
               "model_of_rate_heterogeneity_line2_value","state_frequencies")
    # Make a list of the output rows for the first output dataframe
    op <- c(op1,"amino-acid",op2,op3,num_vals,op4,mrh1,mrh2_name,mrh2,sf1)
    op_df <- data.frame(names,op, stringsAsFactors = FALSE)
    names(op_df) <- c("parameter","value")
    
    # Check whether a gamma matrix is needed
    mrh1_check <- gsub(" ","",mrh1)
    if (mrh1_check=="Uniform"){
      # If the model for rate heterogeneity is uniform, don't need to create a matrix for discrete gamma rate categories
      g_df <- "Uniform"
    } else {
      #Create the matrix for discrete gamma categories
      g_start <- grep(" Category",iq_file)+1 # get the index for the first line of the gamma categories matrix
      empty   <- which(iq_file=="") # get indexes of all empty lines
      empty   <- empty[empty>g_start] # get empty lines above gamma categories matrix
      g_end   <- empty[1]-1 # get end index for gamma categories matrix (one less than next empty line)
      end_line <- iq_file[g_end]
      # if the end isn't an empty line, subtract one from the end count 
      # to exclude lines like "Relative rates are computed as MEAN of the portion of the Gamma distribution falling in the category."
      # to see if this is what's happening, check whether the line starts with a numeric section (i.e. a category for the gamma rate)
      check_line <- length(strsplit(strsplit(end_line, "        " )[[1]][1]," ")[[1]])
      if (check_line > 3){
        # If the check_line is longer than 3 objects, it won't be a group for the gamma categories but an instruction
        # Eg. the relative rates line has 17 objects, because each word in the sentence is an object
        # Instructions can be excluded from the gamma matrix (but categories can't)
        g_end = g_end - 1
      }
      # initialise columns to store data in
      g1 <- c()
      g2 <- c()
      g3 <- c()
      # Iterate through rows in gamma matrix
      for (i in g_start:g_end){
        row <- strsplit(iq_file[[i]],"        ") # split the rows on the long strong of 0's in the middle
        g1 <- c(g1,as.numeric(row[[1]][1])) # add the values to the columns
        g2 <- c(g2,as.numeric(row[[1]][2]))
        g3 <- c(g3,as.numeric(row[[1]][3]))
      }
      g_df <- data.frame(g1,g2,g3, stringsAsFactors = FALSE) # create a dataframe of the information
      names(g_df) <- c("category","relative_rate","proportion") # name the columns 
    }
    
    # Check whether state frequencies are needed
    sf1_squashed <- gsub(" ","",sf1)
    if (sf1_squashed == "(empiricalcountsfromalignment)"){
      # Get starting line for frequencies
      start_ind <- grep("State frequencies:",iq_file) + 2
      # Take the 20 lines containing AA frequencies
      freq_lines <- iq_file[start_ind:(start_ind+19)]
      # Split up the frequency lines into the label and the frequency
      freq_split <- unlist(strsplit(freq_lines,"="))
      # Get the frequency
      freq_nums <- freq_split[c(FALSE,TRUE)]
      # Remove any spaces (from IQTree formatting)
      freq_nums <- gsub(" ","",freq_nums)
      # Get corresponding AA letter
      freq_names <- freq_split[c(TRUE,FALSE)]
      # Remove IQTree formatting
      freq_names <- gsub("pi\\(","",freq_names)
      freq_names <- gsub("\\)","",freq_names)
      freq_names <- gsub(" ","", freq_names)
      # Create a nice dataframe
      f_df <- data.frame("amino_acid" = freq_names,
                         "frequency" = freq_nums,
                         stringsAsFactors = FALSE)
    } else {
      f_df <- "State frequencies from model"
    }
    
    # Create a list of the dataframes - this will be the output
    params <- list(op_df,g_df,f_df)
    # Name the parameters so they're easy to access once you've outputted the data
    names(params) <- c("parameters","gamma_categories", "frequency")
  }
  
  # Now the information has been collected, create an output dataframe
  # make the output dataframe
  return(params)
}


