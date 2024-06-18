# Treelikeness Metrics
#### A comparison of methods for quantifying treelikeness

Caitlin Cherryh

June 2024

***
### Summary
An underlying assumption in phylogenetics is that each site in a loci shares an identical evolutionary history that fits a single bifurcating tree. However, this assumption is broken by biological processes such as introgression or recombination. A number of metrics to exist quantify the deviation from treelikeness, but to our knowledge there is no study comparing the effectiveness of these metrics. To test the metrics, we simulated sequence alignments with two different underlying causes of treelikeness. 

The first approach increases the number of underlying topologies from 1 to 10,000 within a 10KB sequence alignment. When there is 1 underlying tree, the data is perfectly treelike and every site in the alignment shares the same evolutionary history. However, in the most extreme case, there are 10,000 underlying trees, so each site has a different (randomly generated) evolutionary history. 

Our second approach is to simulate a tree and add an introgression event between two taxa. These events can be recent (between two tips that do not form a sister pair) or ancient (occurring between two branches deep in the tree). In this case, the proportion of the alignment varies from 0% (no recombination, the alignment is perfectly treelike) to 50% (half the sites in the alignment have the alternate evolutionary history formed by the introgression event). 

We assess the adequacy of treelikeness metrics to detect changes in treelikeness under these simulation conditions. In addition, we introduce a new treelikeness metric we call the tree proportion, and provide code to calculate the tree proportion for any sequence alignment.  

This github repository contains scripts used to:

1. Simulate data with different underlying levels of treelikeness and incomplete lineage sorting (ILS)
2. Benchmark existing metrics for treelikeness against these simulations
3. Introduce a new metric for treelikeness in phylogenetic datasets, called the tree proportion

If you replicate any part of these analyses or use functions from these scripts, please cite this repository.


#### Contents
+ Scripts
    + All scripts necessary to completely replicate this analysis are included in the `code/` folder
    + Each script includes an overview, a list of necessary parameters or file paths,  and a list of software necessary to run that script
+ Conda enviroment
    + The `environment.yml` file is included to replicate the conda environment used for this project
+ Instructions for replication
    + Instructions for replicating these analyses, along with details about the datasets and software used, are in this `README.md` file.

***
### Instructions to reproduce the analyses:
To fully replicate the analyses, follow these steps:

1. Download and install the software programs and R packages required to run each of the treelikeness metrics
    + IQ-Tree2
    + SplitsTree4 (v4.17.1 or higher)
    + Phylogemetric
    + fast_TIGER
2. Create the conda environment `gene_filtering` using the `environment.yaml` file
3. Prepare simulated alignments
    + The script `code/01_simulations.R` will generate both sets of simulated alignments (random tree simulations and introgression simulations)
    + The random tree simulations mimic decreased treelikeness by increasing the number of evolutionary histories present within an alignment
    + The introgression simulations mimic decreased treelikeness by increasing the proportion of sites within the alignment impacted by a single introgression event
4. Apply treelikeness metrics
    + The script `code/02_apply_treelikeness_metrics.R` will fully replicate our analysis by applying each treelikeness metric to each simulated alignment
    + To apply specific treelikeness metrics to a single alignment, use the functions in the script `code/func_metrics.R`
    + To apply the tree proportion test to a single alignment, use the functions in the script `code/func_tree_proportion.R`
5. Process and analyse results
    + The script `code/03_data_analysis.R` will plot results from the treelikeness metrics for both sets of simulations

***
### Tree proportion

#### Background
The tree proportion metric works by identifying splits within a phylogenetic sequence alignment. A split is a bipartition of the taxa, which indicates there is signal within the alignment (biological or otherwise) supporting separating the set of taxa into two different groups. The presence of incompatible splits means that there is conflicting signal within that data. The tree proportion metric aims to quantify the proportion of information lost in transforming an alignment into a tree by calculating the maximum possible proportion of compatible splits. We find the tree proportion performs well at detecting alignments with multiple distinct topologies, but struggles to detect decreasing treelikeness when alignments contain a single introgression event. However, the tree proportion is a proxy for treelikeness - it does not determine the biological or technical reason for conflict within the alignment.

#### Calculating tree proportion
To calculate the tree proportion, firstly we estimate a phylogenetic network from the input sequence alignment. We then apply a modified version of Kruskal's algorithm to build a greedy tree with the maximum possible sum of split weights. Finally, we find the tree proportion by dividing the sum of split weights in the tree by the sum of split weights in the network.

A tree proportion value of 1 indicates each split in the network is also in the greedy tree. Here, there is no conflicting signal and the dataset is perfectly treelike. Alternatively, a tree proportion of 0.5 indicates that half of the splits in the alignment are conflicting splits. This indicates the level of conflicting signal within the dataset.

#### Trivial splits
There are two ways to calculate the tree proportion: with or without trivial splits. A trivial split is a terminal branch, and terminal branches are present within both the phylogenetic network and the greedy tree. We recommend using the setting `remove_trivial_splits = TRUE` to remove trivial splits, which provides a better picture of the conflicting signal within the alignment of interest.

#### Calculating the tree proportion 
The tree proportion of any alignment can be calculated using the function `tree.proportion`, which is available in the script `code/tree_proportion.R`.

To call the tree proportion function:
```
tree.proportion(alignment_path, sequence_format = "DNA", remove_trivial_splits = TRUE, 
                network.method = "SPLITSTREE" , splitstree_path = NA, dist.ml.model = NA)
```

Function arguments:

- `alignment_path`: valid file path to a sequence alignment (Nexus, Fasta or Phylip formats)
- `sequence_format`: either "DNA" or "AA", depending on your data. Tree proportion has been tested exclusively with DNA data.
- `remove_trivial_splits`: `TRUE` to remove trivial splits (i.e. terminal branches) or FALSE to keep trivial splits. 
    + We recommend the option `remove_trivial_splits = TRUE`
    + All tree proportion values in both sets of simulations were calculated using `remove_trivial_splits = TRUE`
- `network.method`: approach for estimating the phylogenetic network from the provided sequence alignment - either "SPLITSTREE" or "PHANGORN"
    + Options are "SPLITSTREE" (to estimate the phylogenetic network in SplitsTree4) or "PHANGORN" (to estimate the phylogenetic network using the functions `phangorn::dist.ml` and `phangorn::neighborNet`)
    + When `network.method = "SPLITSTREE"`, you must provide the path to the SplitsTree4 java executable (e.g. on my computer this is `/Applications/SplitsTree/SplitsTree.app/Contents/MacOS/JavaApplicationStub`)
    + When `network.method = "PHANGORN"`, a pairwise distance matrix is computed for the input alignment (using `phangorn::dist.ml`), and used as input to estimate a phylogenetic network (using `phangorn::neighbotNet`)
        + You must provide the model used to compute the pairwise distance matrix. Allowed options include "JC69", "F81", and 17 amino acid models
        + See the help page at `?phangorn::dist.ml` for more details about allowable models for this function
        + Extra model parameters (base frequencies, a Q matrix, or the number of rate categories) can be added using the function `calculate.dna.pairwise.distance.matrix()` (located in the `code/tree_proportion.R` script)
        
For example, to apply the tree proportion using SplitsTree4 to estimate the network:
```
tree.proportion("test_alignment.nexus", 
                sequence_format = "DNA", 
                remove_trivial_splits = TRUE, 
                network.method = "SPLITSTREE", 
                splitstree_path = "/Applications/SplitsTree/SplitsTree.app/Contents/MacOS/JavaApplicationStub")
```

Alternatively, to apply the tree proportion using `phangorn::dist.ml` with the F81 model to estimate the network:
```
tree.proportion("test_alignment.nexus", 
                sequence_format = "DNA", 
                remove_trivial_splits = TRUE, 
                network.method = "PHANGORN", 
                dist.ml.model = "F81")
```

#### Network estimation method
We included the "PHANGORN" method of phylogenetic network estimation for users who wanted to apply the tree proportion method without installing and calling a separate software program. However, all simulations for this project used SplitsTree4 to calculate the phylogenetic network. SplitsTree4 is able to estimate phylogenetic networks quickly even for large numbers of taxa (such as our 100 taxa simulations). In our preliminary tests, we found the "PHANGORN" method computationally intractable for large networks, and slower in general than SplitsTree4.

***
### Treelikeness metrics
- **Proportion of explained variance**
    - Note: Referred to as "Cunningham metric" in our manuscript
    - Cunningham, J.P., 1978. "Free trees and bidirectional trees as representations of psychological distance". *Journal of Mathematical Psychology*, 17:165–188. https://doi.org/10.1016/0022-2496(78)90029-9
    - Available in this repository, via the script `code/func_metrics.R` by calling the function `cunningham.test()`
        - Cherryh, C., 2023. "Treelikeness metrics", GitHub repository. https://github.com/caitlinch/treelikeness-metrics (Accessed 6/7/2023)
- **Delta plots**
    - Holland, B.R., Huber, K.T., Dress, A., Moulton, V., 2002. "δ plots: a tool for analyzing phylogenetic distance data". *Molecular Biology and Evolution*, 19:2051–2059. https://doi.org/10.1093/oxfordjournals.molbev.a004030
    - Function `delta.plot()` included in the R package ape
        - Paradis E., Schliep, K., 2019. “ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R.” *Bioinformatics*, 35:526-528. https://doi.org/10.1093/bioinformatics/bty633
- **Likelihood mapping**
    - Strimmer, K., von Haeseler, A., 1997. "Likelihood-mapping: a simple method to visualize phylogenetic content of a sequence alignment". *Proceedings of the National Academy of Science U.S.A*., 94:6815–6819. https://doi.org/10.1073/pnas.94.13.6815
    - Available in IQ-Tree
        - Minh, B.Q., Schmidt, H.A., Chernomor, O., Schrempf, D., Woodhams,  M.D., von Haeseler, A., Lanfear, A., 2020. "IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era". *Molecular Biology and Evolution*, 37:1530-1534.
- **Network treelikeness test**
    - Huson, D.H., Bryant, D., 2006. "Application of phylogenetic networks in evolutionary studies". *Molecular Biology and Evolution*, 23:254-267. https://doi.org/10.1093/molbev/msj030
    - Available in this repository, via the script `code/func_metrics.R` by calling the function `network.treelikeness.test()`. Requires SplitsTree4 (v4.17.1 or above) to run.
        - Cherryh, C., 2023. "Treelikeness metrics", GitHub repository. https://github.com/caitlinch/treelikeness-metrics (Accessed 6/7/2023)
    - SplitsTree4 available from https://uni-tuebingen.de/en/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree/
        - Huson, D.H., Bryant, D., 2006. "Application of phylogenetic networks in evolutionary studies". *Molecular Biology and Evolution*, 23:254-267. https://doi.org/10.1093/molbev/msj030
- **Q-residual**
    - Gray, R.D., Bryant, D., Greenhill, S.J., 2010. "On the shape and fabric of human history". *Philosophical Transactions of the Royal Society B: Biological Sciences*, 365:3923–3933. https://doi.org/10.1098/rstb.2010.0162
    - Available in the Python program Phylogemetric
        - Greenhill, S. J., 2016. "Phylogemetric: A Python library for calculating phylogenetic network metrics". *Journal of Open Source Software*, 1:28. https://doi.org/10.21105/joss.00028
        - Greenhill, S. J., 2016. "phylogemetric", GitHub repository, https://github.com/SimonGreenhill/phylogemetric (Accessed 6/7/2023)
- **Site concordance factors**
    - Minh, B.Q., Hahn, M.W., Lanfear, R., 2020. "New methods to calculate concordance factors for phylogenomic datasets". *Molecular Biology and Evolution*, 9:2727-2733. https://doi.org/10.1093/molbev/msaa106
    - Mo, Y.K., Lanfear, R., Hahn, M.W., Minh, B.Q., 2022. "Updated site concordance factors minimize effects of homoplasy and taxon sampling". *Bioinformatics*, 39:btac741. https://doi.org/10.1093/bioinformatics/btac741
    - Available in IQ-Tree
        - Minh, B.Q., Schmidt, H.A., Chernomor, O., Schrempf, D., Woodhams,  M.D., von Haeseler, A., Lanfear, A., 2020. "IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era". *Molecular Biology and Evolution*, 37:1530-1534. https://doi.org/10.1093/molbev/msaa015 
- **TIGER**
    - Cummins, C.A., McInerney, J.O., 2011. "A Method for Inferring the Rate of Evolution of Homologous Characters that Can Potentially Improve Phylogenetic Inference, Resolve Deep Divergence and Correct Systematic Biases". *Systematic Biology*, 60:833–844. https://doi.org/10.1093/sysbio/syr064
    - Available in the software program fast_TIGER
        - Frandsen, P.B., 2014. "fast-TIGER", Github repository, https://github.com/pbfrandsen/fast_TIGER (Accessed 6/7/2023)

***
### Empirical data
We applied treelikeness metrics to two empirical datasets: 
- **Original dataset:** Metazoa_Choano_RCFV_strict.phy
    - Whelan N.V., Kocot K.M., Moroz T.P., Mukherjee K., Williams P., Paulay G., Moroz L.L., Halanych K.M., 2017. "Ctenophore relationships and their placement as the sister group to all other animals". _Nature Ecology & Evolution_. 1:1737–1746 https://doi.org/10.1038/s41559-017-0331-3
    - Whelan N.V., Kocot K.M., Moroz T.P., Mukherjee K., Williams P., Paulay G., Moroz L.L., Halanych K.M., 2017. "Ctenophora Phylogeny Datasets and Core Orthologs". Dataset. Version 1. https://doi.org/10.6084/m9.figshare.4484138.v1
**Filtered dataset:** 5_Whelan2017MCRS.tar.gz
    - McCarthy, Charley G, P., Mulhair,  P. O., Siu-Ting, K., Creevey, C. J., O’Connell, M. J., 2023. "Improving Orthologous Signal and Model Fit in Datasets Addressing the Root of the Animal Phylogeny", _Molecular Biology and Evolution_, 40:msac276, https://doi.org/10.1093/molbev/msac276
    - McCarthy C. 2022. "ATOLRootStudy", GitHub repository, https://github.com/chmccarthy/ATOLRootStudy (Accessed 18/06/2024)

***
### Citation information
If you replicate any part of these analyses or use functions from these scripts, please cite this repository. Thank you! 

Cherryh, C., 2024. "Treelikeness metrics", GitHub repository. https://github.com/caitlinch/treelikeness-metrics (Accessed 18/6/2024)
