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
* fuzzywuzzy - if it is not installed, try pip install fuzzywuzzy
## Functions:
* it loads your CSV or TXT tables one by one,
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


## How to use:
The merger operates on an interactive way. To run it, just type into your terminal:

python data_cleanup_hausfrau.py

First it asks for your files one by one. Next it will ask for defining the range of the headers. For instance, when your header is only the first row, your start coordinate is 0 and your end coordinate is 1. The first 2 rows of headers would be 0 and 2 , the first 3 rows of headers would be 0 and 3  and so forth.  
Next the field separator is defined. It can be either a tab, or a comma. After that the basic cleanup comes table by table: the header of the table is printed. Here the first row is the index number of the table. This is important for the later editing!  

Then the editing starts: when the column names are the exact same in the given data table, one can choose which columns should be dropped. Also, optimal dropping of unnecessary columns are offered. One can substitute NaNs to zeros as well, like in case of the Prochlorococcus measurements in this subproject.  
Optimal bacterial growth rate can be calculated, just give the index number of the bacterial generation time (Growth rate (per day) = ln2 / generation time).  

Sometimes columns of different tables are slightly different, making the merging impossible. For that, this program performs a fuzzy string matching based on  Levenshtein Distance of a given column name and offers the 3 clostest hits. Then one has to determine which one is the actual match.  
then checks the tables one by one for dupicated and unnecessary rows and offers optional editing. At the end all matching columns will be concatenated to a final table. At the end just for sure the column names are checked again and if the occurennce of a column name is less than the number of input tables, it will be dropped.  
One can give the name of the merge file and the tables are concatenated and inner merged.

If you are not sure how to use it, please check the testrun.txt file!
## Limitations:
 * It is only been tested for Python 2.7
 * Handling categorical data is not solved yet,
 * it cannot handle different table structures. Headers should have the same structure as well,
 * to reduce resetting, it is recommended to have identical headers in every table,
 * it can be problematic to merge and manually set the parameters for hundreds of tables, or really noisiy ones. It is recommended to do the merge with subsets, then merge them together,
 * unique columns with no common match in all input tables will be dropped at the end step
 * it takes what it gets. For instance, when a CSV table is exported from Libreoffice, digits can be lost. Before running this merger, make sure you exported your tables right!
 
## Resources:

https://github.com/seatgeek/fuzzywuzzy

## Have questions or suggestions?
Please contact me at emese.x.szabo@gmail.com!
If you experience bugs or errors, please send me your output message of your terminal, optionally the first 10 lines of your files!
