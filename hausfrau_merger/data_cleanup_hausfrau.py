#!/usr/bin/env/python

# Merging more or less similar CSV format files for later imputation/clustering.

# Details:
# title           :data_cleanup_hausfrau_merge.py
# author          :Emese Xochitl Szabo
# email:	  :emese.szabo@uni-oldenburg.de
# date            :15/02/2021
# version         :0.2
# license         :GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
# usage           :python data_cleanup_hausfrau_merge.py
# notes           :
# python_version  :2.7.16

# Follow the interactive instructions: first a modified header will show your columns, their names, index number,
# and the first 5 rows. Duplicate removal, optional dropping and calculations will be based on the index numbers!

# Every contibution is welcome! To do so please fork the  hausfrau_merger on GitHub, create your chnages locally, and  send a pull request!

import datetime
import operator
import math
import scipy
import numpy as np
import pandas as pd
import os
from itertools import chain
from collections import Counter
import difflib
import fuzzywuzzy
from fuzzywuzzy import fuzz
from fuzzywuzzy import process


now = datetime.datetime.now()

print now

class load_data():

    def __init__(self,input,rows,delim):
        self.loadfile = input
        self.rows = rows
        self.delim=delim
        self.load_table()

    def load_table(self):

        if len(self.rows) > 1:

            self.loadlist = pd.read_csv(self.loadfile, header=self.rows, sep=self.delim)
            self.loadlist.columns = self.loadlist.columns.map('_'.join)
            self.loadlist = self.loadlist.replace(r'^\s*$', np.nan, regex=True)
            self.loadlist = self.loadlist.dropna(how="all", axis=[0,1])
            self.loadlist = self.loadlist.T.drop_duplicates().T
            #self.loadlist = self.loadlist.replace('\W', '', regex=True).convert_objects(convert_numeric=True) # Removes formatting, has to be tested!

            self.loadlist.columns = self.loadlist.columns.str.strip().str.replace('Unnamed: \d+_level_\d+_', '') # when there is no other name, it will retain this. Later at the matching, it will be offered to be merged.
            self.columns = self.loadlist.columns
            self.tableindex = pd.DataFrame(columns=range(len(self.columns)), index=[0])
            self.tableindex.loc[0] = self.loadlist.columns 
            self.tableindex = self.tableindex.astype('str') 

            # for visualization of th indexes and values:
            self.visual = self.loadlist[0:5]
            self.visual.loc[-1] = list(range(len(self.columns)))
            self.visual.index = self.visual.index +1
            self.visual = self.visual.sort_index()

        else:
            self.loadlist = pd.read_csv(self.loadfile, header=rows, sep =delim)
            self.loadlist = self.loadlist.replace(r'^\s*$', np.nan, regex=True)
            self.loadlist = self.loadlist.dropna(how="all", axis=[0,1])
            self.loadlist = self.loadlist.T.drop_duplicates().T
            self.loadlist.columns = self.loadlist.columns.str.strip().str.replace('Unnamed: \d+_level_\d+_', '') # when there is no other name, it will retain this. Later at the matching, it will be offered to be me$
            #self.loadlist = self.loadlist.replace('\W', '', regex=True).convert_objects(convert_numeric=True) # Removes formatting, hast to be tested!

            self.columns = self.loadlist.columns
            self.tableindex = pd.DataFrame(columns=range(len(self.columns)), index=[0])
            self.tableindex.loc[0] = self.loadlist.columns 
            self.tableindex = self.tableindex.astype('str') 

            # for visualization of th indexes and values:
            self.visual = self.loadlist[0:5]
            self.visual.loc[-1] = list(range(len(self.columns)))
            self.visual.index = self.visual.index +1
            self.visual = self.visual.sort_index()

        return self.loadlist, self.columns

inputs=[]  # gonna be raw input

### Welcome message ###

print "Welcome to the Hausfrau merger. This program will help you to clean and concatenate your data tables! Before we start, please consider the following points: \n"

#print recommendations

print "What can this merger do? \n \
- it loads your comma or, tab-separated tables one by one,\n \
- keeps your column names, even with special characters,\n \
- checks column formatting,\n \
- deletes empty rows and columns,\n \
- deletes non-numeric and non-alphabetical characters,\n \
- detects and based on your wish, drops unnecessary duplicated columns,\n \
- additionally, optional column dropping per table is possible too,\n \
- it is able to substitute missing values with zeros,\n \
- it is specialized to Oceanic metadata,therefore growth rate calculation based on bacterial generation time is possible,\n \
- at the merging step, it will detect possibly identical columns and ask you if you would still rename and use them in the merge,\n \
- and finally, it will merge your tables into one final CSV output, based on the columns represented in all your tables\n"

print " What cannot this merger do?\n \
 - It cannot handle different table structures. Headers should have the same structure as well,\n \
 - To reduce resetting, it is recommended to have identical headers in every table,\n \
 - It can be problematic to merge and manually set the parameters for hundreds of tables, or really noisiy ones. It is recommended to do the merge with subsets, then merge them together,\n \
 - Unique columns with no common match in all input tables will be dropped at the end step,\n \
 - It takes what it gets. For instance, when a CSV table is exportedfrom Libreoffice, digits can be lost. Before running this merger, make sure you exported your tables right! \n"

print "If you feel like you made a mistake, no worries! Just break out by hitting ctrl + c!\n"
print "Let's begin!"
print "Please type your input files. Once you are done, hit Enter again!"

sentinel = '' # ends when this string is seen

for line in iter(raw_input, sentinel):
    inputs.append(line)

print inputs

print "Please define the headers of your file. For instance, when your header is only the first row, your start coordinate is 0 and your end coordinate is 1.\n\
The first 2 rows of headers would be 0 and 2 , the first 3 rows of headers would be 0 and 3  and so forth..."
row_start = raw_input("Your header starts at: ")
row_end = raw_input("Your header ends at: ")

print "Please define your field separator (tabulator is \\t, comma is ,):"
delim = raw_input("Your delimiter is: ")

rows = range(int(row_start),int(row_end))
objects = [load_data(input,rows,delim) for input in inputs]
pd.set_option('display.max_rows', 30)
tables = [obj.loadlist for obj in objects]
allcols = []

for obj in objects:
    print "\nProcessing ", obj.loadfile
    cols = [obj.columns]
    cols2 = list(chain.from_iterable(cols))
    pd.set_option('display.max_columns', len(cols2))
    counts = Counter(cols2)
    print "\nHeaders, index numbers of columns and first 5 rows of input: "
    print obj.visual

    droplist=[]
    nanlist=[]
    baclist=[]

    # Detect duplicated columns based on column names.
    for x in counts:
        key = x
        value = counts[key]
        if value > 1:
            print  "\nDuplicated columns were detected in case of ", key, value, "\n"
            answer=raw_input("Would you like to delete duplicated columns? Please answer Yes or No! ").lower()
            if answer=="yes":
                print "The following duplicated columns were found. Which one would you like to keep?"
                print obj.tableindex[~(obj.tableindex != key)].dropna(how="all", axis=[0,1])   # indexes and dup. names in 
                print "Please type the index numbers of the colums to be dropped. Once you are ready, hit Enter!"
                for line in iter(raw_input, sentinel):
                    droplist.append(int(line))
            else: continue
        else: continue

    # Optimal dropping
    answer2=raw_input("\nWould you like to drop any other column? Please answer Yes or No! ").lower()
    if answer2=="yes":
        print "Please type the index numbers of the colums to be dropped. Once you are ready, hit Enter!"
        for line in iter(raw_input, sentinel):
            droplist.append(int(line))
    else: pass

    # Offer NaN to zero:
    answer3=raw_input("\nWould you like to substitute missing column values with zeros, for instance in case of Prochlorococcus measurements? Please answer Yes or No! ").lower()
    if answer3=="yes":
        print "Please type the index numbers of the colums to be substituted. Once you are ready, hit Enter!"
        for line in iter(raw_input, sentinel):
            nanlist.append(int(line))

        updated = obj.loadlist.iloc[:, nanlist].fillna(0)
        obj.loadlist.update(updated)
    else: pass

    # Offer growth rate - instert column at the END
    answer4=raw_input("\nWould you like to calculate growth rate based on bacterial generation time? Please answer Yes or No! ").lower()
    if answer4=="yes":
        print "Please type the index number of the bacterial generation time. Once you are ready, hit Enter!"
        for line in iter(raw_input, sentinel):
            baclist.append(int(line))

        #add new column
        newcol = math.log(2)/obj.loadlist.iloc[:, baclist]
        obj.loadlist['Growth rate'] = newcol

    # Actual dropping
    droplist.sort()
    column_numbers = [x for x in range(obj.loadlist.shape[1]) if x not in droplist] 
    obj.loadlist = obj.loadlist.iloc[:, column_numbers]

    allcols.append(obj.loadlist.columns)

allcols2 = list(chain.from_iterable(allcols))

# test the intersection
allin = set(allcols[0]).intersection(*map(set,allcols))

full_list = set(allcols2)

#unique
diff = {element for element in full_list if element not in allin}
print "\nI have found non-matching column names:\n"
print diff

# String match to save columns:

matchlist = dict()

if len(diff) > 1:
    print "Some headers are not exactly matching. To find possible typos, a similarity analysis is going to be performed. The program tries to find the 3 closest match.\n"
    diff = list(diff)

    for i in diff:

        i_mod = fuzzywuzzy.utils.full_process(i, force_ascii=True)
        closest = process.extract(i_mod,diff,limit=4)
        if len(closest) > 1:
            print "Clostest 3 matches to ", i, ": ", repr(closest[1:]).decode('string-escape') #.encode("utf-8")
            print "Please choose for ", i, "substituting header name (can be one of the clostest matches, or ",i, " as well)"
            newname = raw_input("New header: ")
            matchlist[i] = newname

print "The following headers will be corrected:\n", matchlist

print matchlist
newcols = []
for obj in objects:
    #correction
    obj.loadlist = obj.loadlist.rename(columns=matchlist)
    newcols.append(obj.loadlist.columns)
    #duplicate rows?
    obj.loadlist =  obj.loadlist.drop_duplicates()

newcols2 = list(chain.from_iterable(newcols))
newcounts = Counter(newcols2)
print "If the number of occurrence of header is less than the number of input tables, it will be dropped!"
print newcounts
tables = [obj.loadlist for obj in objects]
mergetable = pd.concat(tables,join='inner',ignore_index=True,sort=False)
print "\nMerged table:"
print mergetable

output=raw_input("Name of the mergefile (csv): ")
try:
    os.remove(output)
except:
    pass

mergetable.to_csv(output, index=False, header= True)
#df.types
print "\nGut gemacht!"
