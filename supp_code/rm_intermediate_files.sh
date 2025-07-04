#!/bin/bash

## Quick script to remove unneeded (intermediate) files generated during simulations


## Supp A1
# Specify exp_A1 directory
a1_dir=/mnt/data/dayhoff/home/u5348329/treelikeness_supp/suppA/exp_A1/
# Delete files `{sample}_output_alignment.fa_Splitstree_output.nex`
find "$a1_dir" -type f -name "*_output_alignment.fa_Splitstree_output.nex" -print
find "$a1_dir" -type f -name "*_output_alignment.fa_Splitstree_output.nex" -delete
# Delete files `{sample}_output_alignment.fa_converted.nex`
find "$a1_dir" -type f -name "*_output_alignment.fa_converted.nex" -print
find "$a1_dir" -type f -name "*_output_alignment.fa_converted.nex" -delete
# Delete files `{sample}_output_alignment.fa.mldist`
find "$a1_dir" -type f -name "*_output_alignment.fa.mldist" -print
find "$a1_dir" -type f -name "*_output_alignment.fa.mldist" -delete
# Delete files `{sample}_output_alignment.fa.ckp.gz`
find "$a1_dir" -type f -name "*_output_alignment.fa.ckp.gz" -print
find "$a1_dir" -type f -name "*_output_alignment.fa.ckp.gz" -delete
# Delete files `{sample}_output_alignment.fa.bionj`
find "$a1_dir" -type f -name "*_output_alignment.fa.bionj" -print
find "$a1_dir" -type f -name "*_output_alignment.fa.bionj" -delete


## Supp A2
# Specify exp_A2 directory
a2_dir=/mnt/data/dayhoff/home/u5348329/treelikeness_supp/suppA/exp_A2
# Delete files `{sample}_output_alignment.fa_Splitstree_output.nex`
find "$a2_dir" -type f -name "*_output_alignment.fa_Splitstree_output.nex" -print
find "$a2_dir" -type f -name "*_output_alignment.fa_Splitstree_output.nex" -delete
# Delete files `{sample}_output_alignment.fa_converted.nex`
find "$a2_dir" -type f -name "*_output_alignment.fa_converted.nex" -print
find "$a2_dir" -type f -name "*_output_alignment.fa_converted.nex" -delete
# Delete files `{sample}_output_alignment.fa.mldist`
find "$a2_dir" -type f -name "*_output_alignment.fa.mldist" -print
find "$a2_dir" -type f -name "*_output_alignment.fa.mldist" -delete
# Delete files `{sample}_output_alignment.fa.ckp.gz`
find "$a2_dir" -type f -name "*_output_alignment.fa.ckp.gz" -print
find "$a2_dir" -type f -name "*_output_alignment.fa.ckp.gz" -delete
# Delete files `{sample}_output_alignment.fa.bionj`
find "$a2_dir" -type f -name "*_output_alignment.fa.bionj" -print
find "$a2_dir" -type f -name "*_output_alignment.fa.bionj" -delete


## Supp B1
# Specify exp_B1 directory
b1_dir=/mnt/data/dayhoff/home/u5348329/treelikeness_supp/suppB/exp_B1/
# Delete files `{sample}_output_alignment.fa.lmap.svg`
find "$b1_dir" -type f -name "*_output_alignment.fa.lmap.svg" -print
find "$b1_dir" -type f -name "*_output_alignment.fa.lmap.svg" -delete
# Delete files `{sample}_output_alignment.fa.lmap.eps`
find "$b1_dir" -type f -name "*_output_alignment.fa.lmap.eps" -print
find "$b1_dir" -type f -name "*_output_alignment.fa.lmap.eps" -delete
# Delete files `{sample}_output_alignment.fa.mldist`
find "$b1_dir" -type f -name "*_output_alignment.fa.mldist" -print
find "$b1_dir" -type f -name "*_output_alignment.fa.mldist" -delete
# Delete files `{sample}_output_alignment.fa.ckp.gz`
find "$b1_dir" -type f -name "*_output_alignment.fa.ckp.gz" -print
find "$b1_dir" -type f -name "*_output_alignment.fa.ckp.gz" -delete
# Delete files `{sample}_output_alignment.fa.bionj`
find "$b1_dir" -type f -name "*_output_alignment.fa.bionj" -print
find "$b1_dir" -type f -name "*_output_alignment.fa.bionj" -delete


## Supp B2
# Specify exp_B2 directory
b2_dir=/mnt/data/dayhoff/home/u5348329/treelikeness_supp/suppB/exp_B2/
# Delete files `{sample}_output_alignment.fa.lmap.svg`
find "$b2_dir" -type f -name "*_output_alignment.fa.lmap.svg" -print
find "$b2_dir" -type f -name "*_output_alignment.fa.lmap.svg" -delete
# Delete files `{sample}_output_alignment.fa.lmap.eps`
find "$b2_dir" -type f -name "*_output_alignment.fa.lmap.eps" -print
find "$b2_dir" -type f -name "*_output_alignment.fa.lmap.eps" -delete
# Delete files `{sample}_output_alignment.fa.mldist`
find "$b2_dir" -type f -name "*_output_alignment.fa.mldist" -print
find "$b2_dir" -type f -name "*_output_alignment.fa.mldist" -delete
# Delete files `{sample}_output_alignment.fa.ckp.gz`
find "$b2_dir" -type f -name "*_output_alignment.fa.ckp.gz" -print
find "$b2_dir" -type f -name "*_output_alignment.fa.ckp.gz" -delete
# Delete files `{sample}_output_alignment.fa.bionj`
find "$b2_dir" -type f -name "*_output_alignment.fa.bionj" -print
find "$b2_dir" -type f -name "*_output_alignment.fa.bionj" -delete


