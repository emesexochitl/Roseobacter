## Intro - What is TDA and how does it work?
The basic idea behind Topological DAta Analysis (TDA) is that data has shape and shape has information (Carlsson, 2009). This comes in handy when we face complex biological data sets with multiple patterns hidden in. Usually in such cases it is really hard to detect these complex, often localized patterns. TDA gives an intuitive insight how these patterns are organized.
All of the classical ordination methods suffer from the disadvantage of pairwise comparisons and necessary exclusion of other less powerful, but important drivers of the ecological processes, or in case of multiple factor analysis, they have to apply constraints (ref and check), further weakening its significance. Nonetheless, because of the interconnectivity, it is hard to pick the best possible hypothesis and how to rank collected information a priori. 

Explain more wih the Ayasdi concept
## Dependencies:

The script is written in Python 3. Because in the time of creating this pipeline many functions were not avaiable in one package, lot of the functions come from various sources. The most important are the following:

* GUDHI: Persistent homology, simplex tree. In order to run the complete script, first the GUDHI environment has to be installed via conda:  
conda install -c conda-forge gudhi  
Then it has to be activated: 
conda activate gudhi  
* numpy, scipy, pandas, Counter for data formatting and special operations
* kmapper for constucting the simplicial complex: https://kepler-mapper.scikit-tda.org/en/latest/
* sklearn and different parts of it to construct lenses
* Mapper: http://danifold.net/mapper/installation/
* sklearn_tda: https://github.com/MathieuCarriere/sklearn-tda
* statmapper: bootstrapping and statistical evaulation. https://github.com/MathieuCarriere/statmapper  

## How to use it: ##
  
As a general step it is recommended first to test how the different lenses perform and turn off the bootstrapping part of the statistical evaulation. As a rule of thumb, a good lens would give a nice, structured inverse image with observable topological elements, such as branches, loops and even separate clusters. If the lens looks nice, one can check the dispersion of values in the coloring step of the graph. It is especially interesting to observe the categorical values represented as piecharts. If everything checks out, one can run the complete script and check the significance of the topological elements.

But what does exacly happening when the script is executed? It begins with loading the environment. Then the input tables are loaded. Here one can choose what type of table they want to use as the base of the graph (in this case it was an OTU/ASV table to see how biodiversity mirrors biogeographical characteristics). The same table can be used for coloring as well, however in this case a separate table was used (environmental parameters). If there are categorical values in this input, they will be transformed to numerical values, but it has to be manually given. After that, one can choose and test the lenses.  
When the simplicial complex is constructed, gain and resolution can be estimated (see statmapper documentation). Based on experience, it gives rather too simple graphs, so as a solution, gain is set to the maximum, 0.4 (personal recommendation of Mathieu Carrière). KMeans clustering with a minimum amount of two centroids is performed in this step.  
After the construction of the simplicial complex, in order to analyze persistent homology it is transformed to a simplex tree, then to persistence diagrams. One type of it is the barcode diagram, where horizontal lines represent the number and dimenson of topological elements in the complex. That can be used to clarify and compare different TDA runs.
networkX

## Functions:
 Lenses:
 * Isolation forest: Custom 1-D lens with Isolation Forest - detecting anomalies.
 * Isomap: 1-D lens with Isomap - to disperse and reduce dimensions.
 * l2norm: Normalization - get some distance with l2norm (check resorurce for better explanation!)
 * MDS: Multidimensional scaling -  display dissimiliarities. It runs a metric MDS.
 * mean: row mean.
 * median: row median.
 * Eccenticity: how close a point lies to the center of the cloud.
 * DTM: Distance to measure - estimating the persistent homology of S directly.
 * L-inf: centrality using Bray-Curtis disssimilarity matrix- assigns to each point the distance to the point most distant from it. Note: many other distance measurements can be used, but BC usually performs the best for ecological samples.
 * PCA: reducing dimensions, finding axes of greatest variance.
 * Laplacian: graph Laplacian of the neighborhood graph.
 * SVD: Singular value decomposition with the first eigenvector - compressing information to smaller space without radical information loss.
 * choosing a column as lens
 * combine two lenses (not applied to the subproject)

## Limitations:
 * if different data sets are used, Min-Max normalization of the data would be recommened for later statistical tests.
 * when using row mean or median, in general lenses don't perform well and a measurement-dependent bias is introduced.
 * It was only tested for the subproject
 
## Resources:
Carlsson, Gunnar. (2009). Topology and Data. Bulletin of The American Mathematical Society - BULL AMER MATH SOC. 46. 255-308. 10.1090/S0273-0979-09-01249-X. 
https://gudhi.inria.fr/  
https://mathieucarriere.github.io/website/publis/publis.html  
https://github.com/MathieuCarriere  
Daniel Müller  
Ayasdi  
