# caitlinch/treelikeness_metrics/code/func_parametric_bootstraps.R
# Caitlin Cherryh 2023

# This file contains functions for generating parametric bootstraps from an empirical alignment



#### Model extraction functions ####
get.simulation.parameters <- function(dotiqtree_file){
  # Given a .iqtree file, this function will extract the relevant parameters
  
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


