#!/usr/bin/env python3

# import necessary modules
import sys
from collections import defaultdict
import pickle

##### This script will be used to build a pickle file that contains a dictionary of the conversion code and the tuple of flags that it corresponds to.

### This section defines the function from the conversion script that was used to convert the flags into a flag code.
def convertFlagsToFlagCode(isReferenceFlag, hasAnchorFlag, lengthFlag, lengthMatchFlag, hasMultipleAnchorsFlag, isQueryFlag, generalErrorFlag):
    flagCode = 0
    if isReferenceFlag == True:
        flagCode += 1
    if hasAnchorFlag == True:
        flagCode += 2
    if lengthFlag == True:
        flagCode += 4
    if lengthMatchFlag == True:
        flagCode += 8
    if hasMultipleAnchorsFlag == True:
        flagCode += 16
    if isQueryFlag == True:
        flagCode += 32
    if generalErrorFlag == True:
        flagCode += 64
    return flagCode

### This section will initialize the dictionary that will be used to store the conversion code and the tuple of flags that it corresponds to.
print("initializing dictionary")
conversionCodeDictionary = defaultdict(list)

### This section will loop through all of the possible flag combinations and add them to the dictionary.
print("building conversion code dictionary")
for isReferenceFlag in [True, False]:
    for hasAnchorFlag in [True, False]:
        for lengthFlag in [True, False]:
            for lengthMatchFlag in [True, False]:
                for hasMultipleAnchorsFlag in [True, False]:
                    for isQueryFlag in [True, False]:
                        for generalErrorFlag in [True, False]:
                            flagCode = convertFlagsToFlagCode(isReferenceFlag, hasAnchorFlag, lengthFlag, lengthMatchFlag, hasMultipleAnchorsFlag, isQueryFlag, generalErrorFlag)
                            flags = (isReferenceFlag, hasAnchorFlag, lengthFlag, lengthMatchFlag, hasMultipleAnchorsFlag, isQueryFlag, generalErrorFlag)
                            conversionCodeDictionary[flagCode] = flags


print("Pickling conversion code dictionary.")
### This section will pickle the dictionary.
pickle.dump(conversionCodeDictionary, open("conversionCodeDictionary.pkl", "wb"))

print("Done.")