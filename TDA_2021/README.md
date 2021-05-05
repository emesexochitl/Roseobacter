## Intro - What is TDA and how does it work?
Intro.

## Dependencies:
In order to run the complete script, first the GUDHI environment has to be installed via the anaconda xxx:
conda install gudhi
Then it has to be activated:
conda activate gudhi


* 
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

## Limitations:
 * if different data sets are used, Min-Max normalization of the data would be recommened for later statistical tests.
 * when using row mean or median, in general lenses don't perform well and a measurement-dependent bias is introduced.
 * It was only tested for the subproject
 
## Resources:
GUDHI
Daniel Müller
Ayasdi
