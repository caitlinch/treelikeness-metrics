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

## Checking progress of NTLT
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
status <- c(
paste0("Exp A1 complete: ", num_A1_complete, " of 3000 (", pc_A1_complete, "%)"),
paste0("Exp A2 complete: ", num_A2_complete, " of 10,000(", pc_A2_complete, "%)")
)
print(status)
```


