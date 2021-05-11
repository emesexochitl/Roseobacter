# Pst-TDA statistical testing.

# Details:
# title           :Stat_test_postTDA.R 
# author          :Emese Xochitl Szabo
# email:          :emese.szabo@uni-oldenburg.de
# date            :01/03/2021
# version         :0.1
# license         :GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
# usage           :Rscript Stat_test_postTDA.R
# notes           :better run in Rstudio and modify to the unique needs (see most commented lines). One can substitute the Wilcoxon-tests to KS-tests.
# r_version       :R.3.6.3

library(dplyr)
library(data.table)
library(tidyverse)

### Environmental parameters ###
rm(list = ls())
#setwd("/home/emese/Desktop/Roseo_R/Roseobacter_2021/") # <- set home folder
tmp <- read.csv("Linf_BC_10_04_new_both_species_scaled_filtered_1lens_stations_node_table.txt", sep=",")

nodelist = c(12,6,15,7) # <- here choose the node, or list of nodes
looplist = colnames(tmp) # go through all variables
print(nodelist)
for (i in colnames(tmp)) {
  if (class(tmp[[i]]) == "integer") {
    selected_data = tmp %>% select(Node, i) %>% dplyr::filter(Node %in% nodelist)
    bg_data = tmp %>% select(Node, i) %>% dplyr::filter(!Node %in% nodelist)
    X <-sprintf("%s,%.3f,%.3f", i, wilcox.test(selected_data[[i]], bg_data[[i]])$statistic, wilcox.test(selected_data[[i]], bg_data[[i]])$p.value)
    if (wilcox.test(selected_data[[i]], bg_data[[i]])$p.value < 0.05) {
      print(X)
    }
    
  }
  else if (class(tmp[[i]]) == "numeric") {
    selected_data = tmp %>% select(Node, i) %>% dplyr::filter(Node %in% nodelist)
    bg_data = tmp %>% select(Node, i) %>% dplyr::filter(!Node %in% nodelist)
    X <- sprintf("%s,%.3f,%.3f", i, wilcox.test(selected_data[[i]], bg_data[[i]])$statistic, wilcox.test(selected_data[[i]], bg_data[[i]])$p.value)
    if (wilcox.test(selected_data[[i]], bg_data[[i]])$p.value < 0.05) {
      print(X)
    }
  }
}

# ######################
# 
# lshap <- lapply(roseo_full_stat, shapiro.test)
# lres <- sapply(lshap, `[`, c("statistic","p.value"))