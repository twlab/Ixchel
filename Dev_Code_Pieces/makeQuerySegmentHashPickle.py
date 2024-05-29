#!/usr/bin/env python3

# import necessary modules
import sys
from collections import defaultdict
import pickle


def rec_dd():
    return defaultdict(rec_dd)

bed_dict = rec_dd()

INPUTFILE = sys.argv[1]
OUTPUTFILE = sys.argv[2]

print("Input file: ")
print(INPUTFILE)


with open(INPUTFILE) as f:
    for line in f:
        L = line.strip().split()
        SEGMENTID = L[1]
        SEQUENCELENGTH = len(L[2])
        STABLESOURCE = ""
        STABLEOFFSET = ""
        #STABLESOURCE = L[3].split(":")[2]
        #STABLEOFFSET = L[4].split(":")[2]


        if bed_dict[SEGMENTID]["SegmentLength"]:
            # Not first entry for segment
            print("ERROR 1")
            print(bed_dict[SEGMENTID]["SegmentLength"])
        else:
            # First postion for read
            bed_dict[SEGMENTID]["StableSource"] = STABLESOURCE
            bed_dict[SEGMENTID]["SegmentLength"] = SEQUENCELENGTH
            bed_dict[SEGMENTID]["StableOffset"] = STABLEOFFSET


# Save nested dictionary as pickle file
print("Saving to: ")
print(OUTPUTFILE)


f = open(OUTPUTFILE, "wb")
pickle.dump(bed_dict, f)
f.close()
