#!/usr/bin/env python3

### This script will read in a GraphMethyl File named GraphMethylFile and a pickle file named ConversionPickleFile.
# It will then read in each line of the GraphMethylFile and lookup the corresponding converted position in the ConversionPickleFile using the graphSegmentID and graphSegmentOffset in the GraphMethylFile.
# It will then write out a new file named MethylCFile with the converted positions.

### GraphMethylFile whose format is:
### 1. This column is the segment ID. Named segmentID
### 2. This column is the offset of the segment. Named segmentOffset
### 3. This column is the sense of the cytosing. Named sense
### 4. This column is the context of the cytosine. Named context
### 5. This column is the number of unmethylated reads. Named unmethylatedReads
### 6. This column is the number of methylated reads. Named methylatedReads
### 7. This column is the total number of reads. Named totalReads
### 8. This column is the fraction of reads that support the methylated fraction of the cytosine. Named methylatedFraction

### ConversionPickleFile whose format is:
### 1. The key is a tuple of the form (graphSegmentID, graphSegmentOffset)
### 2. The value is a tuple of the form (STABLESOURCEID, STABLESOURCESTART, STABLESOURCEEND, FLAGCODE)

### MethylCFile whose format is:
### 1. This column is the converted stable sequence ID. Named STABLESOURCEID
### 2. This column is the starting position of the cytosine in the converted stable sequence. Named STABLESOURCESTART
### 3. This column is the ending position of the cytosine in the converted stable sequence. Named STABLESOURCEEND
### 4. This column is the context of the cytosine. Named CONTEXT
### 5. This column is the fraction of reads that support the methylated fraction of the cytosine. Named METHYLATION
### 6. This column is the sense of the cytosine. Named SENSE
### 7. This column is the total read coverage of the cytosine. Named COVERAGE
### 8. This column is the flag code for the cytosine. Named FLAGCODE


### import necessary modules
import sys
from collections import defaultdict
import pickle

### This section will import the ConversionPickleFile and create a dictionary of the data in the ConversionPickleFile named ConversionDictionary.
ConversionDictionary = defaultdict(list)

print('Reading in ConversionPickleFile')
with open(sys.argv[2], 'rb') as ConversionPickleFile:
    ConversionDictionary = pickle.load(ConversionPickleFile)
ConversionPickleFile.close()

### This section will read in each line of the GraphMethylFile and lookup the corresponding converted position in the ConversionPickleFile using the graphSegmentID and graphSegmentOffset in the GraphMethylFile.
### It will then write out a new file named MethylCFile with the converted positions.
print('Reading in GraphMethylFile and writing out MethylCFile with error codes')
with open(sys.argv[1], 'r') as GraphMethylFile, open(sys.argv[3], 'w', buffering=10000000) as MethylCFile:
    print("Starting to process sites...")
    for line in GraphMethylFile:
        line = line.strip().split('\t')
        conversion = ConversionDictionary[(line[0], line[1])]
        # check if conversion is the right length if it is not then line should be printed with an error code of 128. Then continue to next line.
        if len(conversion) != 4:
            #print("Error: conversion is not the right length. Conversion is: " + str(conversion) + " and line is: " + str(line))
            flagCode = 128
            MethylCFile.write("NA" + '\t' + "NA" + '\t' + "NA" + '\t' + line[3] + '\t' + line[7] + '\t' + line[2] + '\t' + line[6] + '\t' + str(flagCode)+":"+line[0]+":"+line[1] + '\n')
        else:
            flagCode = conversion[3]
            MethylCFile.write(conversion[0] + '\t' + conversion[1] + '\t' + conversion[2] + '\t' + line[3] + '\t' + line[7] + '\t' + line[2] + '\t' + line[6] + '\t' + str(flagCode)+":"+line[0]+":"+line[1] + '\n')
GraphMethylFile.close()



# with open(sys.argv[1], 'r') as GraphMethylFile, open(sys.argv[3], 'w') as MethylCFile:
#     for line in GraphMethylFile:
#         line = line.strip().split('\t')
#         conversion = ConversionDictionary[(line[0], line[1])]
#         flagCode = convertFlagsToFlagCode(conversion[3], conversion[4], conversion[5], conversion[6], conversion[7], conversion[8], conversion[9])
#         MethylCFile.write(conversion[0] + '\t' + conversion[1] + '\t' + conversion[2] + '\t' + line[3] + '\t' + line[7] + '\t' + line[2] + '\t' + line[6] + '\t' + str(flagCode) + '\n')
# GraphMethylFile.close()

