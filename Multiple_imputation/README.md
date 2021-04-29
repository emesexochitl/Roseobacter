## Missing values? Here is how to find and use them!
When one works with complex data tables coming from various resources, mistakes and errors can happen. Sometimes it is manifesting in missing values and in order to save as much valuable measurement as possible, one can use multiple imputation.
## Dependencies:
To perform the multiple imputation, the method ... of missMDA is used, but the complete script depends on more packages:
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

## Functions:
* First, some analytics are performed: normality-check, missing value co-occurrence, ... These measurements will help to better understand the data and help to make good choices later. .. are visulaized as well
* Sesonf, multiple imputation is performed on nonscaled data, then on scaled data
* Min-Max normlization
* In the last part of the R script, several ordination methods are included, as an alternative to analyze ecologivcal data sets
* 
* 
* 
* 
* 
* 
* 

## Limitations:
 * It is only been tested for R 3.6.3
 * Usually imputing missing values more than 20 % will lead less trustworthy results. Handle it with precaution!
 * It was only tested for the subproject, where different data tables were imputed one by one, then merged to one table (see the Hausfrau folder)
 * 
## Resources:
Josse J, Husson F (2016). “missMDA: A Package for Handling Missing Values in Multivariate Data Analysis.” Journal of Statistical Software, 70(1), 1–31. doi: 10.18637/jss.
