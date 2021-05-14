## Intro - What is TDA and how does it work?
The basic idea behind Topological Data Analysis (TDA) is that data has shape and shape has information (Carlsson, 2009). This comes in handy when we face complex biological data sets with multiple patterns hidden in. Usually in such cases it is really hard to detect these complex, often localized patterns. All classical ordination methods suffer from the disadvantage of pairwise comparisons and necessary exclusion of other less powerful, but important drivers of the ecological processes, or in case of multiple factor analysis, they have to apply constraints (ref and check), further weakening its significance. Nonetheless, because of the interconnectivity, it is hard to pick the best possible hypothesis and how to rank collected information a priori.  TDA gives an intuitive insight how these patterns are organized.

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

For the detailed depencencies, please check the script!  

## How to use it: ##
  
As a general step it is recommended first to test how the different lenses perform and turn off the bootstrapping part of the statistical evaulation. As a rule of thumb, a good lens would give a nice, structured inverse image with observable topological elements, such as branches, loops and even separate clusters. If the lens looks nice, one can check the dispersion of values in the coloring step of the graph. It is especially interesting to observe the categorical values represented as piecharts. If everything checks out, one can run the complete script and check the significance of the topological elements.  

But what does exacly happening when the script is executed? It begins with loading the environment. Then the input tables are loaded. Here one can choose what type of table they want to use as the base of the graph (in this case it was an OTU/ASV table to see how biodiversity mirrors biogeographical characteristics). The same table can be used for coloring as well, however in this case a separate table was used (environmental parameters). If there are categorical values in this input, they will be transformed to numerical values, but it has to be manually given. After that, one can choose and test the lenses.  

When the simplicial complex is constructed, gain and resolution can be estimated (see statmapper documentation). Based on experience, it gives rather too simple graphs, so as a solution, gain is set to the maximum, 0.4 (personal recommendation of Mathieu Carrière). KMeans clustering with a minimum amount of two centroids is performed in this step.  
After the construction of the simplicial complex, in order to analyze persistent homology it is transformed to a simplex tree, then to persistence diagrams. One type of it is the barcode diagram, where horizontal lines represent the number and dimenson of topological elements in the complex. That can be used to clarify and compare different TDA runs.

<p float="left">
<img src="https://user-images.githubusercontent.com/14163953/117679982-0966ec80-b1b1-11eb-8b13-2993308d2aca.jpg" width="500" height="auto" />
<img src="https://user-images.githubusercontent.com/14163953/117679987-09ff8300-b1b1-11eb-81d0-f45eac010ce2.jpg"  width="500" height="auto" />
</p>

Other analytic steps are performed as well, but before that, the graph is transformed into a NetworkX graph. This allows to do comparative analysis by "freezing" the structure. However, a classic KeplerMapper HTML output is also generated, together with a PDF version.  

<p float="left">
<img src="https://user-images.githubusercontent.com/14163953/117680785-c78a7600-b1b1-11eb-87c4-1dd3695f2932.png" width="500" height="auto" />
<img src="https://user-images.githubusercontent.com/14163953/117680791-c8230c80-b1b1-11eb-9c56-32711ae3607c.jpg"  width="500" height="auto" />
</p>

In the next step, the statistical evaulation takes place. It is recommened to have 100 times bootstrapping with 95% CI. As a result, the significant topological elements (connected components, loops, and up- or downbranches) will be colored yellow. Beware, when it tries to evaulate the down/upbranch it analyzes the values IN the graph. After that, two grey colored plots are generated as well, to show either the node numbers, or names. 

<p float="left">
<img src="https://user-images.githubusercontent.com/14163953/117681623-952d4880-b1b2-11eb-91a7-63bbcd4e9ea2.jpg" width="300" height="300" />
<img src="https://user-images.githubusercontent.com/14163953/117681611-92325800-b1b2-11eb-9df1-2230bf0ee45f.jpg" width="300" height="300" />
<img src="https://user-images.githubusercontent.com/14163953/117681598-8fcffe00-b1b2-11eb-8170-f0cec0333f19.jpg"  width="300" height="300" />
</p>

Sometimes when graphs get more complex, other visualizations are useful as well. Betwenness and Eigenvector-centrality are plotted for this reason. With betweenness we can see the shortest way, or the level of connections of the nodes, with the Eigenvector version we see the influence of a node has on the network.

<p float="left">
<img src="https://user-images.githubusercontent.com/14163953/117682028-f5bc8580-b1b2-11eb-88fc-9253cf8e8ea4.jpg" width="500" height="auto" />
<img src="https://user-images.githubusercontent.com/14163953/117682019-f523ef00-b1b2-11eb-844c-c4e8f64c97d0.jpg"  width="500" height="auto" />
</p>

In the last part, coloring of the graph takes place. In case of TDA not only the shape, but the value representation has important meaning. While shapes can mean cyclic nature (loops), or anomalies (branches), coloration can highlight other hidden subparts of the graph. Since two input files were given, two types of coloration happens. One colors based on the average value of a given environmental parameter per node. Exception is the categorical value, where modal values are used. In bot cases, node size is normalized for better, balaced visualization, but the values displayed as colors and barplots are not normalized. For better visualization, categorical values are plotted as piecharts to see the exact distribition of values in the nodes. Important here is to give manually the column names of the categorical values! The outputs are saved in PDF files and the continous value plots are compiled into a video. 
  
<p float="left">
<img src="https://user-images.githubusercontent.com/14163953/117807783-5eac0800-b25c-11eb-9b49-f647db9b4f27.png" width="300" height="auto" />
<img src="https://user-images.githubusercontent.com/14163953/117807781-5e137180-b25c-11eb-904c-649f4e889d2b.png"  width="300" height="auto" />
<img src="https://user-images.githubusercontent.com/14163953/117807774-5ce24480-b25c-11eb-92a3-a8dbb9cdf50f.png"  width="300" height="auto" />
</p>

The same procedure is repeated for the OTU/ASV table, but numerical identifiers are added to the IDs because of the Unknown/same names. Node values are extracted for later statistical tests. 

## Functions:
 Lenses:
 * Isolation forest: Custom 1-D lens with Isolation Forest - detecting anomalies.
 * Isomap: 1-D lens with Isomap - to disperse and reduce dimensions.
 * l2norm: Normalization - get some distance with l2norm (check resorurce for better explanation!)
 * MDS: Multidimensional scaling -  display dissimiliarities. It runs a metric MDS.
 * mean: row mean.
 * median: row median.
 * Eccentricity: how close a point lies to the center of the cloud.
 * DTM: Distance to measure - estimating the persistent homology of S directly.
 * L-inf: centrality using Bray-Curtis disssimilarity matrix- assigns to each point the distance to the point most distant from it. Note: many other distance measurements can be used, but BC usually performs the best for ecological samples.
 * PCA: reducing dimensions, finding axes of greatest variance.
 * Laplacian: graph Laplacian of the neighborhood graph.
 * SVD: Singular value decomposition with the first eigenvector - compressing information to smaller space without radical information loss.
 * choosing a column as lens
 * combine two lenses (not applied to the subproject)

## Limitations:
 * If different data sets are used, Min-Max normalization of the data would be recommened for later statistical tests.
 * With statmapper, only relatively simple graphs can be generated.
 * When using row mean or median, in general lenses don't perform well and a measurement-dependent bias is introduced.
 * I mostly focused on the connected components and loops at the significance test. Hopefully in the future there will be tools for more sophisticated analysis.
 * It was only tested for the subproject
 
## Resources:
Carlsson, Gunnar. (2009). Topology and Data. Bulletin of The American Mathematical Society - BULL AMER MATH SOC. 46. 255-308. 10.1090/S0273-0979-09-01249-X.  
https://gudhi.inria.fr/  
https://mathieucarriere.github.io/website/publis/publis.html  
https://github.com/MathieuCarriere  
http://danifold.net/mapper/index.html   
https://www.ayasdi.com/social/TDAintroduction/  

## Have questions or suggestions?
Please contact me at emese.szabo@uni-oldenburg.de!
If you experience bugs or errors, please send me your output message of your terminal, optionally the first 10 lines of your files!
