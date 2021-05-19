## How to run the TINA/PINA distance measurement script?

As an iteresting way of analyzing oceanograpic data and biodiversity, one can use the functions.community_similarity.R script by Schmith et. al, 2016. The problem is that thanks to an update of the R package ape, the fast.cophenetic.phylo function will fail.  
To avoid it, one has to follow the next steps:  
1. git clone the https://github.com/rec3141/Schmidt_et_al_2016_community_similarity
2. Go to https://github.com/rec3141/Schmidt_et_al_2016_community_similarity/commit/de26dafdd30c4229749c93054a0891227dfca77e?branch=de26dafdd30c4229749c93054a0891227dfca77e&diff=split and change the 	node.ages = ape_node_depth_edge_length(Ntip = Ntip(x), Nnode = x$Nnode, edge = x$edge, Nedge = nrow(x$edge)[1], edge.length = x$edge.length) to node.ages = node.depth.edgelength(x) in your local  functions.community_similarity.R script
3. In the command line,start R or Rstudio and  run the  functions.community_similarity.R:  source(" functions.community_similarity.R", echo=TRUE)
4. Then edit the  prepare.community_similarity.R script: comment out phyloseq at the beginning, and run the script

Because it needs a lot of resurces, a SLURM script and the modified R scipts are also provided in this folder.
