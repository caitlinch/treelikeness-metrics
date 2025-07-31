#!/bin/bash

## Quick script to remove unneeded (intermediate) files generated during simulations


## Exp 1
# Specify exp_1 directory
exp_dir=/mnt/data/dayhoff/home/u5348329/treelikeness_metrics/exp_1/
# Delete files `{sample}_output_alignment.fa_Splitstree_output.nex`
find "$exp_dir" -type f -name "*_output_alignment.fa_Splitstree_output.nex" -print
find "$exp_dir" -type f -name "*_output_alignment.fa_Splitstree_output.nex" -delete
# Delete files `{sample}_output_alignment.fa_Splitstree_NeighborNet_splits.nex`
find "$exp_dir" -type f -name "*_output_alignment.fa_Splitstree_NeighborNet_splits.nex" -print
find "$exp_dir" -type f -name "*_output_alignment.fa_Splitstree_NeighborNet_splits.nex" -delete
# Delete files `{sample}_output_alignment.fa_converted.nex`
find "$exp_dir" -type f -name "*_output_alignment.fa_converted.nex" -print
find "$exp_dir" -type f -name "*_output_alignment.fa_converted.nex" -delete
# Delete files `{sample}_output_alignment.fa_confidence.nexus`
find "$exp_dir" -type f -name "*_output_alignment.fa_confidence.nexus" -print
find "$exp_dir" -type f -name "*_output_alignment.fa_confidence.nexus" -delete
# Delete files `{sample}_output_alignment.fa.ckp.gz`
find "$exp_dir" -type f -name "*_output_alignment.fa.ckp.gz" -print
find "$exp_dir" -type f -name "*_output_alignment.fa.ckp.gz" -delete
# Delete files `{sample}_output_alignment.fa.nex`
find "$exp_dir" -type f -name "*_output_alignment.fa.nex" -print
find "$exp_dir" -type f -name "*_output_alignment.fa.nex" -delete
# Delete files `{sample}_output_alignment.fa.mldist`
find "$exp_dir" -type f -name "*_output_alignment.fa.mldist" -print
find "$exp_dir" -type f -name "*_output_alignment.fa.mldist" -delete
# Delete files `{sample}_output_alignment.fa.bionj`
find "$exp_dir" -type f -name "*_output_alignment.fa.bionj" -print
find "$exp_dir" -type f -name "*_output_alignment.fa.bionj" -delete
# Delete files `{sample}_output_alignment.fa.lmap*`
find "$exp_dir" -type f -name "*_output_alignment.fa.lmap*" -print
find "$exp_dir" -type f -name "*_output_alignment.fa.lmap*" -delete
# Delete files `{sample}_output_alignment.fa.log`
find "$exp_dir" -type f -name "*_output_alignment.fa.log" -print
find "$exp_dir" -type f -name "*_output_alignment.fa.log" -delete
# Delete files `{sample}_output_alignment.fa.iqtree`
find "$exp_dir" -type f -name "*_output_alignment.fa.iqtree" -print
find "$exp_dir" -type f -name "*_output_alignment.fa.iqtree" -delete
# Delete files `{sample}_output_alignment.fa_ParsimonyInformativeSites_only*`
find "$exp_dir" -type f -name "*_output_alignment.fa_ParsimonyInformativeSites_only*" -print
find "$exp_dir" -type f -name "*_output_alignment.fa_ParsimonyInformativeSites_only*" -delete
# Delete files `scfl.*`
find "$exp_dir" -type f -name "scfl.*" -print
find "$exp_dir" -type f -name "scfl.*" -delete



