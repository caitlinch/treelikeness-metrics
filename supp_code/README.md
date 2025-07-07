---
title: "Running supplementary analyses"
author: "Caitlin Cherryh"
date: "2025-06-30"
output: html_document
---

## Activate conda environment:

From home directory on Dayhoff:
```
cd treelikeness_supp
module load mamba
source activate /mnt/data/dayhoff/home/u5348329/treelikeness_supp/envs
```

## Adding packages to environment:
```
module load mamba
mamba activate /mnt/data/dayhoff/home/u5348329/treelikeness_supp/envs
mamba install <package_name>
```

## R packages location (Dayhoff)
`/tmp/RtmpJNryzw/downloaded_packages`

## Checking progress of Supp A: NTLT
In BASH:
```
expA1_dir=/mnt/data/dayhoff/home/u5348329/treelikeness_supp/suppA/exp_A1/
expA1_complete=$(find "$expA1_dir" -name "*NTLT_results.csv" | wc -l)
expA2_dir=/mnt/data/dayhoff/home/u5348329/treelikeness_supp/suppA/exp_A2/
expA2_complete=$(find "$expA2_dir" -name "*NTLT_results.csv" | wc -l)
echo "Exp A1 complete: $expA1_complete of 3000"
echo "Exp A2 complete: $expA2_complete of 10,000"
```

In R:
```
results_dir = "/mnt/data/dayhoff/home/u5348329/treelikeness_supp/suppA/"
expA1_dir <- paste0(results_dir, "exp_A1/")
expA2_dir <- paste0(results_dir, "exp_A2/")
num_A1_complete <- length(grep("NTLT_results.csv", list.files(expA1_dir, recursive = TRUE)))
pc_A1_complete <- round(num_A1_complete/3000*100, digits = 2)
num_A2_complete <- length(grep("NTLT_results.csv", list.files(expA2_dir, recursive = TRUE)))
pc_A2_complete <- round(num_A2_complete/10000*100, digits = 2)
status <- paste0(
"Exp A1 complete: ", 
num_A1_complete, 
" of 3000 (", 
pc_A1_complete, 
"%) ; Exp A2 complete: ", 
num_A2_complete, 
" of 10,000 (", 
pc_A2_complete, 
"%)"
)
print(status)
```

## Checking progress of Supp B: LM

In R:
```{r}
results_dir = "/mnt/data/dayhoff/home/u5348329/treelikeness_supp/suppB/"
expB1_dir <- paste0(results_dir, "exp_B1/")
expB2_dir <- paste0(results_dir, "exp_B2/")
num_B1_complete <- length(grep("LM_results.csv", list.files(expB1_dir, recursive = TRUE)))
pc_B1_complete <- round(num_B1_complete/3750*100, digits = 2)
num_B2_complete <- length(grep("LM_results.csv", list.files(expB2_dir, recursive = TRUE)))
pc_B2_complete <- round(num_B2_complete/5940*100, digits = 2)
status <- paste0(
"Exp B1 complete: ",
num_B1_complete,
" of 3750 (",
pc_B1_complete,
"%) ; Exp B2 complete: ",
num_B2_complete,
" of 5940 (",
pc_B2_complete,
"%)"
)
print(status)
```

## Quick function to check if row output csv files are complete
```
read.complete.csv <- function(file_name){
  if (file.exists(file_name) == TRUE){
    if (file.info(file_name)[["size"]] != 0){
      file_contents <- read.csv(file_name)
      return(file_contents)
    }
  }
}

expA2_ntlt_rows <- lapply(
expA2_ntlt_csvs,
read.complete.csv
)

write.csv(
expA2_ntlt_df,
file = paste0(results_directory, "expA2_NTLT_results_collated.csv"),
row.names = FALSE
)
```


