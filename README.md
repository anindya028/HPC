# HPC
Hierarchical Phylogeny Construction

This program constructs phylogeny in multiple levels from a sequence alignment in fasta format. First, it creates the phylogeny for the moderately similar isolates and reports it as "best_high_resolution_tree.nwk". The highly similar or closely related isolates are reported in a text file named as "<fasta_filename>_groups.txt", where each line contains a comma separated list of the names of closely related isolates in that group.

The phylogeny of each group of closely related isolates is reported as RAxML_bestTree.<group_id>_alignment_sampled for RAxML, <group_id>_alignment_sampled.fasta.treefile for IQtree and <group_id>_FastTree.nwk for FastTree.  

HPC can be called from command line as the following command:

./HPC alignment.fasta -t 16

where "-t 16" means 16 cores are going to be used for computation. 

By default, it will use RAxML for building trees for moderately similar and highly similar (closely related) isolates. If we want to use IQ-TREE or FastTree, we need to use the "-p" option with the first letter of the program name:

./HPC alignment.fasta -t 16 -p F 

This command will result in phylogenies constructed by FastTree for both moderately and highly similar isolates. 
