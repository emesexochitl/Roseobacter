# Roseobacter
## Merge data tables easily!
Welcome to the Hausfrau merger. This program was created to merge and concatenate slightly or more dissimar data tables for further processing.

## Dependencies:
The merger depends on the following Python packages:
import math
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
* deletes empty rows and columns,
* deletes non-numeric and non-alphabetical characters,
* detects and based on your wish, drops unnecessary duplicated columns,
* additionally, optional column dropping per table is possible too,
* it is able to substitute missing values with zeros,
* it is specialized to Oceanic metadata, therefore growth rate calculation based on bacterial generation time is possible,
* at the merging step, it will detect possibly identical columns and ask you if you would still rename and use them in the merge,
- and finally, it will merge your tables into one final CSV output, based on the columns represented in all your tables.
## Limitations:
 * It cannot handle different table structures. Headers should have the same structure as well,
 * To reduce resetting, it is recommended to have identical headers in every table,
 * It can be problematic to merge and manually set the parameters for hundreds of tables, or really noisiy ones. It is recommended to do the merge with subsets, then merge them together,
 * Unique columns with no common match in all input tables will be dropped at the end step
 * It takes what it gets. For instance, when a CSV table is exported from Libreoffice, digits can be lost. Before running this merger, make sure you exported your tables right!
