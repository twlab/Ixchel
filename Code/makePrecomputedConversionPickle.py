#!/usr/bin/env python3

# import necessary modules
import sys
from collections import defaultdict
import pickle

### This script will read in a file PreComputedPositionsFile. It will take the flags and calculate a conversionCode for each. I will then create a dictionary of the data in the PreComputedPositionsFile named ConversionDictionary. It will then pickle ConversionDictionary and save it as a pickle file with the name in PickleFile.

### PreComputedPositionsFile whose format is:
### 1. This column is the converted stable sequence ID. Named STABLESOURCEID
### 2. This column is the starting position of the cytosine in the converted stable sequence. Named STABLESOURCESTART
### 3. This column is the ending position of the cytosine in the converted stable sequence. Named STABLESOURCEEND
### 4. This column is the context of the cytosine, but this column isn't needed for this script. Named CONTEXT
### 5. This is the number of reads that support the methylated fraction of the cytosine, but this column isn't needed for this script. Named METHYLATION
### 6. This is the sense of the cytosine, but this column isn't needed for this script (I think). Named SENSE
### 7. This is the total read coverage of the cytosine, but this column isn't needed for this script. Named COVERAGE
### 8. This is the flagCode of the cytosine. Named FLAGCODE
### 15. This is the graph segment ID. Named graphSegmentID
### 16. This is the offset of the graph segment. Named graphSegmentOffset

### ConversionDictionary whose format is:
### 1. The key is a tuple of the form (graphSegmentID, graphSegmentOffset)
### 2. The value is a tuple of the form (STABLESOURCEID, STABLESOURCESTART, STABLESOURCEEND, isReferenceFlag, hasAnchorFlag, lengthFlag, lengthMatchFlag, hasMultipleAnchorsFlag, isQueryFlag, generalErrorFlag)


### This section will import the PreComputedPositionsFile and create a dictionary of the data in the PreComputedPositionsFile named ConversionDictionary. It will print a message when it has processed 1% increments of the file.
ConversionDictionary = defaultdict(list)
print('Reading in PreComputedPositionsFile')
with open(sys.argv[1], 'r') as PreComputedPositionsFile:
    for line in PreComputedPositionsFile:
        line = line.strip().split('\t')
        ConversionDictionary[(line[8], line[9])] = (line[0], line[1], line[2], line[7])
        if len(ConversionDictionary) % 1000000 == 0:
            print('Processed ' + str(len(ConversionDictionary)) + ' lines so far.')

PreComputedPositionsFile.close()

print('Saving ConversionDictionary as a pickle file')
### This section will pickle ConversionDictionary and save it as a pickle file with the name in PickleFile.
with open(sys.argv[2], 'wb') as PickleFile:
    pickle.dump(ConversionDictionary, PickleFile, protocol=pickle.HIGHEST_PROTOCOL)
PickleFile.close()



