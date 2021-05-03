## Missing values? Here is how to find and use them!
When one works with complex data tables coming from various resources, mistakes and errors can happen. Sometimes it is manifesting in missing values and in order to save as much valuable measurement as possible, one can use multiple imputation.
## Dependencies:
To perform the multiple imputation, the method of performing principal component methods on incomplete data to estimate parameters of missMDA is used, but the complete script depends on more packages:
* ggplot2
* lattice
* reshape2
* naniar
* missMDA
* FactoMineR
* WIM
* factoextra
* corrplot
* data.table
* heatmaply
* vegan
* MASS
* extrafont
* dplyr
* tibble
* divo
* gtools
* vegan3d
* tidyverse

## Functions:
 First, the table is loaded and formatted (check hidden characters, select and reformat samples)
 Analytics are performed: normality-check (Shapiro-test), missing value co-occurrence, ... These measurements will help to better understand the data and help to make good choices later. These steps are visualized as well:
 
* Normality plot: value distributions.

![normality_plot](https://user-images.githubusercontent.com/14163953/116690671-03eff200-a9ba-11eb-9b45-f129ad2070bc.png)

* Matrix plot: visualizing missing values co-occurrences. Missing values are marked as red, and can indicate systematic errors.

![matrixplot](https://user-images.githubusercontent.com/14163953/116692086-1a974880-a9bc-11eb-87c0-7dfdbbb4d4cf.png)


Second, multiple imputation is performed on scaled data for visualization purposes. This includes estimating the  number of dimensions used in the reconstruction formula, then generate the imputed data sets with the MIPCA function using the number of dimensions previously calculated, finally plot the results:

* Combination plot: pattern visualization of missing and imputed values.

![combinations](https://user-images.githubusercontent.com/14163953/116693905-d48fb400-a9be-11eb-997c-60a0b8f8da29.png)

* Correlation matrix plots: visualize positive pairwise correlations of values (negative correlations are not trustworthy!). They come in "upper", "hclust" and "mixed" flavors.

<p float="left">
<img src="https://user-images.githubusercontent.com/14163953/116694281-613a7200-a9bf-11eb-955f-7f38aa90b40c.png" width="300" height="auto" />
<img src="https://user-images.githubusercontent.com/14163953/116694286-61d30880-a9bf-11eb-8bda-67d989a33bea.png"  width="300" height="auto" />
<img src="https://user-images.githubusercontent.com/14163953/116694288-626b9f00-a9bf-11eb-9c1c-1b4411655f1d.png"  width="300" height="auto" />
</p>

* Heatmap: interactive correlation matrix plot with dendrogram.

![interctive_corrplot](https://user-images.githubusercontent.com/14163953/116696945-d3608600-a9c2-11eb-8db6-2e5e164fbed3.png)

* Multiple imputation plots: they show  the procrustean representation of the individuals, dimensions of the PCA, projection of the individuals as supplementary individuals and the projection of the variables as supplementary variables, in this order. Pay attention to the amount of imputed values and inertia on the plots!



After that, multiple imputation is repeated on nonscale data, for one type of ooutput. However, it is important to note for the scaled data imputation that Min-Max normlization is performed. It is based on the fact that most datasets are not normally distributed(see previous Shapiro-test) and this way of normalization provides upper and lower boundaries.

In the last part of the R script, several ordination methods are included, as an alternative to analyze ecological data sets

* Visualizations
 

## Limitations:
 * It is only been tested for R 3.6.3
 * Usually imputing missing values more than 20 % will lead less trustworthy results. Handle it with precaution!
 * It was only tested for the subproject, where different data tables were imputed one by one, then merged to one table (see the Hausfrau folder)
 * 
## Resources:
Josse J, Husson F (2016). “missMDA: A Package for Handling Missing Values in Multivariate Data Analysis.” Journal of Statistical Software, 70(1), 1–31. doi: 10.18637/jss.
