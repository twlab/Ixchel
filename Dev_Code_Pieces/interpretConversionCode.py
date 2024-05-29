#!/usr/bin/env python3

# import necessary modules
import sys
from collections import defaultdict
import pickle

##### This script will be used to interpret the conversion code that is outputted from the conversion script.
# It will do this by loading the conversion code dictionary pickle file and then looking up the conversion code in the dictionary.

### This section reads in the code to be converted from the command line.
conversionCode = int(sys.argv[1])

### This section will load the conversion code dictionary pickle file.
conversionCodeDictionary = pickle.load(open("conversionCodeDictionary.pkl", "rb"))

### This section will look up the conversion code in the dictionary and print the corresponding tuple of flags.
flags = conversionCodeDictionary[conversionCode]

### This section will convert the tuple of flags into a tab delimited string.
flagString = "\t".join([str(flag) for flag in flags])

### This section will print the tab delimited string.
print(flagString)

