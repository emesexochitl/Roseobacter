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
* VIM
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
 First, the table is loaded and formatted (check hidden characters, select and reformat samples), then analytics are performed: normality-check (Shapiro-test) and  missing value co-occurrence.These measurements will help to better understand the data and help to make good choices later. These steps are visualized as well:
 
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

* Multiple imputation plots: they show  the procrustean representation of the individuals, dimensions of the PCA and the projection of the variables as supplementary variables, in this order. Pay attention to the amount of imputed values and inertia on the plots!

<p float="left">
<img src="https://user-images.githubusercontent.com/14163953/117017372-d5448500-acf3-11eb-8721-722d27b3816e.jpg" width="300" height="auto" />
<img src="https://user-images.githubusercontent.com/14163953/117017367-d4abee80-acf3-11eb-9f9d-772675c7c360.jpg"  width="300" height="auto" />
<img src="https://user-images.githubusercontent.com/14163953/117017352-d1b0fe00-acf3-11eb-8d33-b546390cdc9a.jpg"  width="300" height="auto" />
</p>

After that, multiple imputation is repeated on nonscale data, for one type of output. However, it is important to note for TDA Min-Max normalization was performed on imputed, nonscaled data. The main reason was to be able to compare different samplings of origin. Choosing the Min-Max normalization comes in handy when datasets are not normally distributed (see previous Shapiro-test) and this way of normalization provides upper and lower boundaries.

In the last part of the R script, several ordination methods are included, as an alternative to analyze ecological data sets. They weren't tested for this project.

## Limitations:
 * It is only been tested for R 3.6.3
 * Usually imputing missing values more than 20 % will lead less trustworthy results. Handle it with precaution!
 * It was only tested for the subproject, where different data tables were imputed one by one, Min-Max normalized, then merged to one table (see the Hausfrau folder)
 
## Resources:
Josse J, Husson F (2016). “missMDA: A Package for Handling Missing Values in Multivariate Data Analysis.” Journal of Statistical Software, 70(1), 1–31. doi: 10.18637/jss.

https://datasharkie.com/how-to-normalize-data-in-r/

https://www.codecademy.com/articles/normalization

## Have questions or suggestions?
Please contact me at emese.szabo@uni-oldenburg.de!
If you experience bugs or errors, please send me your output message of your terminal, optionally the first 10 lines of your files!
