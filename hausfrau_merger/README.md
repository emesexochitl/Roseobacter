# Hausfrau merger README
## Merge data tables easily!
Welcome to the Hausfrau merger! This program was created to merge and concatenate slightly or more dissimar data tables for further processing.

## Dependencies:
The merger depends on the following Python packages:
* math
* scipy
* numpy
* pandas
* os
* itertools
* collections
* difflib
* fuzzywuzzy 
## Functions:
* it loads your CSV tables one by one,
* keeps your column names, even with special characters,
* checks column formatting,
* deletes whitespaces,
* deletes empty rows and columns,
* deletes non-numeric and special characters,
* detects and based on your wish, drops unnecessary duplicated columns,
* additionally, optional column dropping per table is possible too,
* it is able to substitute missing values with zeros,
* it is specialized to Oceanic metadata, therefore growth rate calculation based on bacterial generation time is possible,
* at the merging step, it will detect possibly identical columns and ask you if you would still rename and use them in the merge,
* and finally, it will merge your tables into one final CSV output, based on the columns represented in all your tables,
* you can choose the name of the output file.
## Limitations:
 * Handling categorical data is not solved yet,
 * it cannot handle different table structures. Headers should have the same structure as well,
 * to reduce resetting, it is recommended to have identical headers in every table,
 * it can be problematic to merge and manually set the parameters for hundreds of tables, or really noisiy ones. It is recommended to do the merge with subsets, then merge them together,
 * unique columns with no common match in all input tables will be dropped at the end step
 * it takes what it gets. For instance, when a CSV table is exported from Libreoffice, digits can be lost. Before running this merger, make sure you exported your tables right!
 ## How to use:
The merger operates on an interactive way. To run it, just type into your terminal:

python data_cleanup_hausfrau.py

First it asks for your files one by one. Next it will ask for defining your headers, then checks the tables one by one for dupicated and unnecessary rows and offers optional editing. At the end all matching columns will be concatenated to a final table.
If you are not sure how to use it, please check the testrun.txt file!
