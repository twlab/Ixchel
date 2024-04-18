#!/usr/bin/env python3

# import necessary modules
import sys
from collections import defaultdict
import pickle


def rec_dd():
    return defaultdict(rec_dd)


bed_dict = rec_dd()
doubleanchor_dict = rec_dd()

INPUTFILE = sys.argv[1]
OUTPUTFILE = sys.argv[2]
DOUBLEANCHORFILE = sys.argv[3]

print("Input file: ")
print(INPUTFILE)


with open(INPUTFILE) as f:
    for line in f:
        L = line.strip().split()
        REFSEGMENTID = L[1]
        QUERYSEGMENTID = L[3]
        REFSTRAND = L[2]
        QUERYSTRAND = L[4]
        TAG = L[5]

        # It looks like for the purposes of merging pickles we have to make the outer most layer ALIGNMENT
        if bed_dict[QUERYSEGMENTID]["RefSegmentID"]:
            # Not first anchor for downstream segment. This is a double anchor. This needs to be fixed, but I am not sure how to do it.
            #print(["Double Anchor detected", QUERYSEGMENTID])
            doubleanchor_dict[REFSEGMENTID]["QuerySegmentID"] = QUERYSEGMENTID
            bed_dict[QUERYSEGMENTID]["RefSegmentID"] = REFSEGMENTID
            bed_dict[QUERYSEGMENTID]["RefStrand"] = REFSTRAND
            bed_dict[QUERYSEGMENTID]["QueryStrand"] = QUERYSTRAND
            bed_dict[QUERYSEGMENTID]["Tag"] = TAG
        else:
            # First position for read
            bed_dict[QUERYSEGMENTID]["RefSegmentID"] = REFSEGMENTID
            bed_dict[QUERYSEGMENTID]["RefStrand"] = REFSTRAND
            bed_dict[QUERYSEGMENTID]["QueryStrand"] = QUERYSTRAND
            bed_dict[QUERYSEGMENTID]["Tag"] = TAG


#print(bed_dict["37372"]["RefSegmentID"])
#print(bed_dict["37372"]["RefStrand"])
#print(bed_dict["37372"]["QueryStrand"])
#print(bed_dict["37372"]["Tag"])


# Save nested dictionary as pickle file
print("Saving to: ")
print(OUTPUTFILE)

f = open(OUTPUTFILE, "wb")
pickle.dump(bed_dict, f)
f.close()

print(DOUBLEANCHORFILE)

f = open(DOUBLEANCHORFILE, "wb")
pickle.dump(doubleanchor_dict, f)
f.close()


