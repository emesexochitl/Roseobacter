
## How to test TDA topological elements for significance?

In order to see if the given topological element, which can be a loop, branch or even one node is significant, here is a small R script to test it.  
## Dependencies:
The merger depends on the following R packages:  
* tidyverse
* dplyr
* data.table

## Functions:
* It needs the node table output from the TDA python script as an input: * \_node_table.csv
* Then to extract values, just type the list of node numbers (it has to be inspected visually based on the TDA graph, see TDA_2021!)
* The list of chosen nodes then will be tested against the leftover data (background), variable by variable
* It will print out the statistically significant variables' (p-val < 0.05) names, test statistic and p-value
* The script will automatically check and convert categorical values 


## Limitations:
 * To test with an other statistical test, one has to exchange the test in the script manually (easy in RStudio)
 * Statistical test is two sided
 
 ## How to use:

## Have questions or suggestions?
Please contact me at emese.szabo@uni-oldenburg.de!
If you experience bugs or errors, please send me your output message of your terminal, optionally the first 10 lines of your files!

