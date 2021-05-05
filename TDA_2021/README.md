## Intro - What is TDA and how does it work?
The basic idea behind Topological DAta Analysis (TDA) is that data has shape and shape has informtion (Carlsson, 2009). This comes in handy when we face complex biological data sets with multiple patterns hidden in. Usually insuch cases it is really hard to detect these complex, often localized patterns. TDA gives an intuitive insight how these patterns are organized

## Dependencies:
In order to run the complete script, first the GUDHI environment has to be installed via the anaconda xxx:  
conda install gudhi  
Then it has to be activated:  
conda activate gudhi  

## How to use it: ##

Begin with loading the environment, then edit the scipt according to your needs. It is recommended first to test how the different lenses perform and turn off the bootstrapping part of the stratistical evaulation. As a rule of thumb, a good lens would give a nice, structured inverse image, with observable topological elements, such as branches, loops and even separate clusters. If the lens looks nice, one can check the dispersion of  values in the coloring step of the graph. It is especially interesting to observe the categorical values represented as piecharts.

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
GUDHI
Daniel Müller
Ayasdi
