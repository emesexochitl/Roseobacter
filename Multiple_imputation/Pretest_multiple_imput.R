# Performing multiple imputation to save data.

# Details:
# title           :Pretest_multiple_imput.R 
# author          :Emese Xochitl Szabo
# email:	        :emese.szabo@uni-oldenburg.de
# date            :26/02/2021
# version         :0.1
# license         :GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
# usage           :Rscript data_cleanup_hausfrau_merge.py
# notes           :better run in Rstudio and modify to the unique needs (see most commented lines)
# r_version       :R.3.6.3


library("ggplot2")
library("lattice")
library("reshape2")
library("naniar")
library(missMDA)
library(FactoMineR)
library("VIM")
library("factoextra")
library("corrplot")
library(data.table)
library(heatmaply)
library(vegan)
library(MASS)
library(extrafont)
loadfonts()
library(dplyr)
library(tibble)
library(gtools)
library(vegan3d)
library(divo)
library(tidyverse)

#setwd("") <- working library
pdf("Meinhard_eco_pacific_scaled.pdf")  
roseo_full <- read.csv("Meta_Data_Pacific_Comparison.txt", header = T,sep="\t",as.is = T, fileEncoding="UTF-8", check.names = F)
roseo_full <- roseo_full %>% filter(str_detect(Sample_ID, "DV", negate = TRUE)) # for the special Deep Vienna samples
roseo_full$Size_Fraction <- as.numeric(roseo_full$Size_Fraction)

## Normality distribution check and visualization ##

roseo_full_stat <-roseo_full[,-c(1,2,7,36)]
lshap <- lapply(roseo_full_stat, shapiro.test) # KS
lres <- sapply(lshap, `[`, c("statistic","p.value"))
lres

roseo_full_scale <-scale(roseo_full_stat)
df <- as.data.frame(roseo_full_scale)
df.m <- melt(df)
ggplot(df.m) + geom_freqpoly(aes(x = value,y = ..density.., colour = variable)) +  theme(legend.position="none", panel.background = element_blank())
matrixplot(roseo_full_scale, cex.axis=.25) 

## Starting multiple imputationon scaled data based on the missMDA method (visualization)##
nb <- estim_ncpPCA(roseo_full_scale,method.cv = "Kfold", verbose = FALSE)
res.comp <- imputePCA(roseo_full_scale, ncp = nb$ncp)
res.pca <- PCA(res.comp$completeObs, ncp = nb$ncp, graph=FALSE)
cor.mat <- round(cor(res.comp$completeObs),2)
summary(res.pca,nbelements=Inf, file="PCA_summary_full_scale.txt")
res<-summary(aggr(roseo_full_scale, sortVar=TRUE))$combinations
corrplot(cor.mat, type="upper", order="hclust",tl.col="black", tl.srt=45, tl.cex=0.6, title="Correlation matrix of filtered data", mar=c(0,0,1,0)) # ?
corrplot(cor.mat, order="hclust",tl.col="black", tl.srt=45, tl.cex=0.6, title="Correlation matrix of filtered data", mar=c(0,0,1,0), addrect = 7) #?
corrplot.mixed(cor.mat, number.cex=0.4, tl.pos="lt", tl.col="black", tl.cex=0.7) # ?
my_cor <- cor(cor.mat)
heatmaply_cor(my_cor, scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "red", high = "blue"), margins = c(40, 40))

## Visualizing multiple imputation ##
resMI <- MIPCA(roseo_full_scale,ncp=nb$ncp) # ncp was already calculated!
#plot(resMI) ### when generating output, this has to be turned off, otherwise it will fail
plot(resMI, choice ="ind.proc", new.plot=FALSE)
plot(resMI, choice ="dim", new.plot=FALSE)
plot(resMI, choice ="ind.supp", new.plot=FALSE)

## Resetting parameters for the last plot. ##
opar <- par()
par(cex = 0.5)
plot(resMI, choice ="var", new.plot=FALSE)
par(opar)

## Multiple imputation with no scaling (one type of output) ## 
roseo_full_nonscale_clean <-roseo_full[,-c(1,2,7,36)]
nb <- estim_ncpPCA(roseo_full_nonscale_clean,method.cv = "Kfold", verbose = FALSE, scale = F)
res.comp <- imputePCA(roseo_full_nonscale_clean, ncp = nb$ncp,scale = F)
export_nonscale_full <- cbind(roseo_full[,c(1,2,7,36)],res.comp$completeObs)
write.table(export_nonscale_full, file = "pacific_nonscale_imputed_meinhard_full.txt", quote = F, sep = "\t")
write.table(res.comp$completeObs, file = "pacific_nonscale_imputed_meinhard_onlyimputed.txt", quote = F, sep ="\t")

## When different measurement series are used, min-max scaling is applied for better coloring and comparative stats of the TDA: (second option for outpput; we used this, run separately for Atlantic and Pacific sets, then merged them) ##
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

# Assemble final output: multiple imputed, originally nonscaled, but Min-Max normalized #
nonscale_table <- as.data.frame(res.comp$completeObs)
scale_table <- apply(nonscale_table,2,normalize)
export_scale_full <- cbind(roseo_full[,c(1,2,7,36)],scale_table)
write.table(export_scale_full, file = "pacific_scale_imputed_meinhard_full.txt", quote = F, sep = "\t")

#############################################################################################################################
#############################################################################################################################

## This part was not used, however can be really handy for later ordination analysis:##

# use the imputed, nonscaled data, then pca, then, NMDS, RDA, PCoA

roseo_vegan <-res.comp$completeObs # nonscaled multiple imputation
#roseo_vegan <- read.table("nonscale_imputed_meinhard_onlyimputed.txt")
roseo_veganspec <- read.csv("TaxCountTab.csv", sep=",") # OTU/ASV tables

##Get the taxonomic level, here the Order ##
roseo_veganspec <- roseo_veganspec[,-c(1,2,3,5,6)]

## Remove NAs ## 
roseo_veganspec$Order <- as.character(roseo_veganspec$Order)
roseo_veganspec<- na.replace(roseo_veganspec, 'Unknown')
na_count <-sapply(roseo_veganspec, function(y) sum(length(which(is.na(y)))))
na_count

## First with DNA ##
roseo_veganspec_dna <- dplyr::select(roseo_veganspec, starts_with("d"))
roseo_veganspec_dna <- cbind(roseo_veganspec[1], roseo_veganspec_dna)
roseo_veganspec_dna <- setDT(roseo_veganspec_dna)[, lapply(.SD, sum), by = Order, .SDcols = d.PS79.178.20m_p.fastq:d.PS79.330.60m_p.fastq]
roseo_veganspec_dna <- column_to_rownames(roseo_veganspec_dna, var="Order")
roseo_veganspec_dna <- roseo_veganspec_dna[ rowSums(roseo_veganspec_dna)!=0, ]
colnames(roseo_veganspec_dna) <- gsub("d.PS79.", "", colnames(roseo_veganspec_dna))
colnames(roseo_veganspec_dna) <- gsub("m_p.fastq", "", colnames(roseo_veganspec_dna))

collist <- colnames(roseo_veganspec_dna)
roseo_vegan_dna <- roseo_vegan[collist,]

## Check if some samples are lost: ##
testlist<-rownames(subset(roseo_vegan,row.names(roseo_vegan) %in% collist))
setdiff(collist,testlist)

roseo_veganspec_dna <- t(roseo_veganspec_dna)

## Check if the subsetting was okay ##
nrow(roseo_vegan_dna)
nrow(roseo_veganspec_dna)

## Rarefaction ## 
dna_specaccum <-specaccum(roseo_veganspec_dna)
plot(dna_specaccum)

div <- diversity(roseo_veganspec_dna, index = "shannon") # add to data, later for TDA!
write.table(div, file="Roseobacter_DNA_alphadiv.txt", quote = F, sep = "\t")

rankindex(scale(roseo_veganspec_dna), roseo_veganspec_dna, c("euc","man","bray","jac","kul"))

vare.dist <- vegdist(roseo_veganspec_dna)
dna_dist <- dist(roseo_vegan)
write.table(as.matrix(dna_dist), file="Roseobacter_euc_distance_params.txt", quote = F, sep = "\t")

# PCA #
vare.pca <- rda(roseo_vegan_dna, scale = TRUE)
plot(vare.pca, scaling = 3)

vare.mds <- metaMDS(roseo_veganspec_dna, trace = TRUE)
saveRDS(vare.mds, "roseo_spec_all_metaMDS_dna.rds")
stressplot(vare.mds, vare.dist)
ordiplot(vare.mds, type = "p", display = "species") 
vare.mds$stress
ordiplot(vare.mds, type = "p", display = "sites")

gof <- goodness(vare.mds)
max.gof <- max(gof)
point.size <- 5/max.gof
sit.sc <- scores(vare.mds)
plot(vare.mds, type="n", main="Goodness of fit")
points(sit.sc,pch=21, cex=gof*point.size)
text(sit.sc, row.names(roseo_veganspec_dna), pos=3, cex=0.7)

## CCA ##
vare.ca <- cca(roseo_veganspec_dna) #, scale = TRUE)
vare.ca
plot(vare.ca, scaling = 3) # scaling here is what to take as an average
chisq.test(roseo_veganspec_dna/sum(roseo_veganspec_dna)) # p-val is not important here!
## Detrented CCA - not recommended ##
vare.dca <- decorana(roseo_veganspec_dna)
vare.dca
plot(vare.dca, display="species") 

## Environmental fitting ##

# have to modify the parameter tables - delete rows, and the species table  - modifying names

ef <- envfit(vare.mds, roseo_vegan_dna, permu = 999)
ef
plot(vare.mds, display = "sites")
plot(ef, p.max = 0.001)
ef2 <- envfit(vare.mds ~ pot.Temperature_C + Bacterial_generation_time_days, roseo_vegan_dna)
plot(vare.mds, display = "sites")
plot(ef2)
tmp <- with(roseo_vegan_dna, ordisurf(vare.mds, pot.Temperature_C, add = TRUE))
with(roseo_vegan_dna, ordisurf(vare.mds, Bacterial_generation_time_days, add = TRUE, col = "green4"))

## create CCA, then envfit ##
ef3 <- envfit(vare.ca, roseo_vegan_dna, permutations = 999, na.rm = T)
plot(vare.ca, display = "sites")
plot(ef3, p.max = 0.001)

## Reformat data tables ##
roseo_veganspec_dna <- as.data.frame(roseo_veganspec_dna)
roseo_vegan_dna <- as.data.frame(roseo_vegan_dna)

vare.cca <- cca(roseo_veganspec_dna ~ pot.Temperature_C + Salinity_psu + Bacterial_generation_time_days, roseo_vegan_dna)
vare.cca
plot(vare.cca)
ef4 <- envfit(vare.cca, roseo_vegan_dna, permutations = 999, na.rm = T)
ef4
plot(vare.cca)#, display = "sites")
plot(ef4)
with(roseo_vegan_dna, ordisurf(vare.cca, pot.Temperature_C, add = TRUE))
with(roseo_vegan_dna, ordisurf(vare.cca, Salinity_psu, add = TRUE, col = "green4"))
with(roseo_vegan_dna, ordisurf(vare.cca, Bacterial_generation_time_days, add = TRUE, col = "blue"))

# with 3D environmental representation: #
ordiplot3d(vare.cca, type = "h")
ordirgl(vare.cca, display = "species", type = "t")
ordirgl(vare.cca, display = "sites", type = "t")

environmental fitting for 3D representation #
ef5 <- envfit(vare.cca, roseo_vegan_dna, permutations = 999, na.rm = T, choices = c(1,2,3))
ordirgl(vare.cca, display = "species", type = "t", envfit=ef5)

## Help building models and choosing the most important environmental factors ##
vare.mod0 <- cca(roseo_veganspec_dna ~ 1,roseo_vegan_dna)
vare.mod1 <- cca(roseo_veganspec_dna ~ .,roseo_vegan_dna)
vare.mod <- step(vare.mod0, scope = formula(vare.mod1), test = "perm")
vare.mod

ef6 <- envfit(vare.mod ~ `POC [microg/l]` + `Amino acid turnover time_[days]`, roseo_vegan_dna) #NMDS is not good for ordisurf since it is ranked! use CA 
plot(vare.mod, display = "sites", scaling = 3)
plot(ef6)
with(roseo_vegan_dna, ordisurf(vare.mod, `POC [microg/l]`, add = TRUE))
with(roseo_vegan_dna, ordisurf(vare.mod, `Amino acid turnover time_[days]`, add = TRUE, col = "green4"))
with(roseo_vegan_dna, ordisurf(vare.mod, Bacterial_generation_time_days, add = TRUE, col = "blue"))
with(roseo_vegan_dna, ordisurf(vare.mod, Eukaryotes_rel_um, add = TRUE, col = "black"))

## Calibration of variables - check if it is curvy, it cannot be used ##
pred <- calibrate(vare.mod)
with(roseo_vegan_dna, plot(Growth_rate_own, pred[,"Growth_rate_own"] - Growth_rate_own, ylab="Prediction Error"))
abline(h=0, col="grey")
plot(mod, display = c("bp", "wa", "lc"))
plot(vare.mod, display = c("bp", "wa", "lc"))
ef7 <- with(roseo_vegan_dna, ordisurf(vare.mod, Eukaryotes_rel_um, display = "lc", add = TRUE))

## 3D represantation ## 
ordirgl(vare.mod, display = "species", type = "t")
ordirgl(vare.mod, display = "sites", type = "t")

## repeat with NDMS, even it is not recommended, because this is the original all modell with abundances! ##
plot(vare.mds, display = "sites")
plot(ef6)
with(roseo_vegan_dna, ordisurf(vare.mds, `POC [microg/l]`, add = TRUE))
with(roseo_vegan_dna, ordisurf(vare.mds, `Amino acid turnover time_[days]`, add = TRUE, col = "green4"))
with(roseo_vegan_dna, ordisurf(vare.mds, Bacterial_generation_time_days, add = TRUE, col = "blue"))
with(roseo_vegan_dna, ordisurf(vare.mds, Eukaryotes_rel_um, add = TRUE, col = "black"))

## Beta diversity ##
betad <- betadiver(roseo_veganspec_dna, "z")
write.table(as.matrix(betad), file="Roseobacter_DNA_betadiv.txt", quote = F, sep = "\t")

adonis(betad ~ `POC [microg/l]`*`Amino acid turnover time_[days]`, roseo_vegan_dna, perm=200)
beta.mod <- with(roseo_vegan_dna, betadisper(betad))
beta.mod <- with(roseo_vegan_dna, betadisper(betad, pot.Temperature_C)) # need category!
plot(beta.mod)
boxplot(beta.mod)
anova(beta.mod)

## Perform Mantel ## 
pc <- prcomp(roseo_vegan_dna, scale = TRUE)
pc<- scores(pc, display = "sites", choices = 1:4)
edis <- vegdist(pc, method = "euclid")
vare.dis <- vegdist(wisconsin(sqrt(roseo_veganspec_dna)))
mantel(vare.dis, edis)
plot(vare.dis, edis)
pc <- scores(pc, choices = 1:2)
pro <- protest(vare.mds, pc) # check here the right mds
plot(pro)

## Classification ##
dis <- vegdist(roseo_veganspec_dna)
clua <- hclust(dis, "average") # middle solution between single lineage and complete
range(dis) # Y axis
plot(clua)
rect.hclust(clua, 5)
grp <- cutree(clua, 5)

## Testing classification ##
boxplot(`POC [microg/l]` ~ grp, data=roseo_vegan_dna, notch = FALSE)
ord <- cca(roseo_vegan_dna)
plot(ord, display = "sites")
ordihull(ord, grp, lty = 2, col = "red")
plot(ord, display="sites"
ordicluster(ord, clua, col="blue")
mst <- spantree(dis, toolong = 1)
plot(mst, ord=ord, pch=21, col = "red", bg = "yellow", type = "t")

## Reorder dendrogram based on first axis of CA ##
wa <- scores(ord, display = "sites", choices = 1)
den <- as.dendrogram(clua)
oden <- reorder(den, wa, mean)
op <- par(mfrow=c(2,1), mar=c(3,5,1,2)+.1)
plot(den)
plot(oden)
par(op)

## Clustering ## 
vegemite(roseo_veganspec_dna, use = oden, zero = "-", scale = "Hill")
tabasco(decostand(roseo_veganspec_dna,"log"), ord)
plotree <- hclust(vegdist(roseo_veganspec_dna), "average")
tabasco(roseo_veganspec_dna, plotree)

dev.off()
dev.off()

## This part can be useful for the environmental dataset ##

## Check which index separates the community the best ##
rankindex(scale(roseo_vegan), roseo_vegan, c("euc","man","bray","jac","kul"))
ggplot(roseo_full, aes(x=roseo_full$Lon, y=roseo_full$`Depth [m]`, colour= roseo_full$`Depth [m]`)) + geom_point() + theme_classic() + scale_y_reverse() + scale_color_gradient2(midpoint = 4000, high = scales::muted("green"), low = scales::muted("blue"), mid = scales::muted("yellow"))
