# HPC
Hierarchical Phylogeny Construction

This program constructs phylogeny in multiple levels from a sequence alignment in fasta format. First, it creates the phylogeny for the moderately similar isolates and reports it as "best_high_resolution_tree.nwk". The closely related isolates are reported in a text file named as "<fast_filename>_groups.txt", where each line contains a comma separated list of the names of closely related isolates in that group.

The phylogeny of each group of closely related isolates is reported as RAxML_bestTree.group_id_alignment_sampled for RAxML, group_id_alignment_sampled.fasta.treefile for IQtree and group_id_FastTree.nwk for FastTree.  

