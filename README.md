# No evidence for functional convergence in animal terrestrialization

The repository contains the data and analysis scripts accompanying the manuscript Kasalo, N., Domazet-Lošo, M., Domazet-Lošo, T. No evidence for functional convergence in animal terrestrialization.

All code was written in RStudio version 2025.05.0+496, R version 4.3.2, on macOS Monterey version 12.0.1.

The repository contains the following files:
1. analyses code.R - the code necessary to perform the analyses presented in the manuscript
2. names.txt - mapping of species names from Wei. et al. Nature 649, 638–646 (2026) to taxIDs
3. nodes.txt - all child-parent pairs in the phylogeny
4. name_mapping.xlsx - updated names for taxa od interest
5. tree_calib.txt - time-calibrated phylogeny from Wei. et al. Nature 649, 638–646 (2026)
6. tree_v2.xml - uncalibrated phylogeny generated based on the calibrated tree from Wei. et al. Nature 649, 638–646 (2026)
7. ps_mapping/gain_c011 - HG - node of origin mapping for Orthogroups from Wei. et al. Nature 649, 638–646 (2026)
8. Phylo_v1_allc_output - numbers of HGs gained at each node in each lineage
9. results - output of permutation analyses (the first row represents the values derived from the 11 terrestrialization nodes, while the other 10,000 rows represent random permutations of nodes):
   a. permutations_noenrichment_gain - shared number of GO terms for each group of nodes (Fig. 1a and Supp. Fig. 1a in our manuscript)
   b. permutations_noenrichment_weighted_gain - node-based average frequency of shared GO terms (Fig. 1b and Supp. Fig. 1b in our manuscript)
   c. permutations_noenrichment_weighted_cluster - cluster-based average frequency of shared GO terms (Fig. 1c and Supp. Fig. 1c in our manuscript)
   d. permutations_gainrate - the gain rate results (Fig. 2 in our manuscript)
    

Other necessary files - these files are too large and are available through FigShare: https://doi.org/10.6084/m9.figshare.33835246
1. Clustering/results_0_11/db_clu_all.tsv - cluster - member mapping, derived from the Orthogroups.txt file from Wei. et al. Nature 649, 638–646 (2026)
2. nature_all_annotated_default.emapper - eggNOG-mapper output for all genomes
3. go.obo - description of GO terms from the Gene Ontology website
