#!/usr/bin/env python3

# import necessary modules
import sys
from collections import defaultdict
import pickle
#import numpy as np



GraphMethylFile,ReferenceSegmentsPickle,QuerySegmentsPickle,LinksPickle,OutputFile,RefOnlyParam,UpstreamOutputFile,DownstreamOutputFile,DoubleAnchorFile = sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8],sys.argv[9]


#### Initialize hashes and files
print("Initialize hashes...")
def rec_dd():
    return defaultdict(rec_dd)

ref_dict = {}
query_dict = {}
link_dict = {}

print("Adding reference segments...")
with open(ReferenceSegmentsPickle, 'rb') as f:
    ref_dict = pickle.load(f)
f.close()

print("Adding query segments...")
with open(QuerySegmentsPickle, 'rb') as f:
    query_dict = pickle.load(f)
f.close()

print("Adding links segments...")
with open(LinksPickle, 'rb') as f:
    link_dict = pickle.load(f)
f.close()

print("Adding double anchor segments...")
with open(DoubleAnchorFile, 'rb') as f:
    doubleanchor_dict = pickle.load(f)
f.close()


print("Input file: ")
print(GraphMethylFile)


### Define function to pull reference segment coordinates
# This setting returns a different format than the other setting!!!
def pullRefOnlyCoords(line):
    L = line.strip().split()
    SEGMENTID = L[0]
    SEGMENTOFFSET = int(L[1])
    SENSE = L[2]
    CONTEXT = L[3]
    # UNMETHYLATED = int(L[4])
    # METHYLATED = int(L[5])
    COVERAGE = int(L[6])
    METHYLATEDFRACTION = L[7]
    if ref_dict[SEGMENTID]:

        STABLESOURCE = ref_dict[SEGMENTID]["StableSource"]
        START = int(ref_dict[SEGMENTID]["StableOffset"]) + SEGMENTOFFSET
        STOP = START + 1
        REFCHECK = True
        return([STABLESOURCE,START,STOP,CONTEXT,METHYLATEDFRACTION,SENSE,COVERAGE,REFCHECK])
    else:
        print("Query Segment")
        return ([False, False, False, False, False, False, False, False])



### This section will define the function that will convert the flags into a flag code.
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



print("Adding upstream links array...")
with open(UpstreamOutputFile, 'rb') as f:
    upstreamkeyarray = pickle.load(f)
f.close()

print("Adding downstream links array...")
with open(DownstreamOutputFile, 'rb') as f:
    downstreamkeyarray = pickle.load(f)
f.close()


def pullAllCoords(line):
    L = line.strip().split()
    SEGMENTID = L[0]
    SEGMENTOFFSET = int(L[1])
    SENSE = L[2]
    CONTEXT = L[3]
    #UNMETHYLATED = int(L[4])
    #METHYLATED = int(L[5])
    COVERAGE = int(L[6])
    METHYLATEDFRACTION = L[7]

    # Test if segment is a reference segment
    REFTEST = ref_dict[SEGMENTID]
    QUERYTEST = query_dict[SEGMENTID]
    FLAG = False
    QUERYCHECK = bool(QUERYTEST)

    # Test if segment is an alt segment and if it is less than 2bp in  length
    if not REFTEST and QUERYTEST:

        try:
            LENGTHTEST = query_dict[SEGMENTID]["SegmentLength"] <= 2
        except TypeError or ValueError:
            LENGTHTEST = False
            print("Error code 1")
            print(SEGMENTID)
            print(query_dict[SEGMENTID]["SegmentLength"])
            print(type(query_dict[SEGMENTID]["SegmentLength"]))
            FLAG = True

    elif REFTEST:
        LENGTHTEST = ref_dict[SEGMENTID]["SegmentLength"] <= 2
        DOUBLEANCHORTEST = False
    else:
        print("Error code 2")
        print("Segment is neith reference nor alt")
        print(SEGMENTID)
        DOUBLEANCHORTEST = False
        LENGTHTEST = False
        FLAG = True

    # Test if segment has an anchor, which is a segment preceding the current one that is a reference segment
    ANCHORTEST = link_dict[SEGMENTID]

    # Test if syntenic reference segment is of a different len than the query segment
    if ANCHORTEST and not REFTEST and QUERYTEST:
        ANCHORSEGMENTID = link_dict[SEGMENTID]["RefSegmentID"]
        SEGMENTSENSE = link_dict[SEGMENTID]["QueryStrand"]
        if doubleanchor_dict[ANCHORSEGMENTID]["QuerySegmentID"]:
            #This tests whether the anchor segment is a double anchor. For some reason the test isn't working.
            SYNTENICLENGHTTEST = False
            FLAG = True
            DOUBLEANCHORTEST = True
        else:
            DOUBLEANCHORTEST = False
            # This should pull the reference segment that is downstream of the anchor segment
            try:
                SyntenicReferenceSegmentID = downstreamkeyarray[upstreamkeyarray.index(ANCHORSEGMENTID)]
                SYNTENICLENGHTTEST = ref_dict[SyntenicReferenceSegmentID]["SegmentLength"] == query_dict[SEGMENTID]["SegmentLength"]
            except ValueError:
                print("Probable double anchor")
                print(ANCHORSEGMENTID)
                SYNTENICLENGHTTEST = False
                FLAG = True
    elif ANCHORTEST and REFTEST:
        SYNTENICLENGHTTEST = False
        DOUBLEANCHORTEST = False
        SEGMENTSENSE = link_dict[SEGMENTID]["RefStrand"]
    elif not ANCHORTEST and REFTEST:
        SYNTENICLENGHTTEST = False
        DOUBLEANCHORTEST = False
        SEGMENTSENSE = "+"
    else:
        SYNTENICLENGHTTEST = False
        DOUBLEANCHORTEST = False
        FLAG = True
        SEGMENTSENSE = "NA"

    # given the set of tests, perform the appropriate lift over for the entry and flag it's status
    if REFTEST:
        # Simplest case, segment is a reference segment
        STABLESOURCE = ref_dict[SEGMENTID]["StableSource"]
        START = int(ref_dict[SEGMENTID]["StableOffset"]) + SEGMENTOFFSET
        # does the sense needs to be considered here?
        STOP = START + 1
        REFCHECK = True
        ANCHORCHECK = bool(ANCHORTEST)
        LENGTHCHECK = LENGTHTEST
        SYNTENICLENGTHCHECK = SYNTENICLENGHTTEST

    elif not REFTEST and ANCHORTEST and LENGTHTEST and SYNTENICLENGHTTEST:
        # Segment is a query segment with an anchor and is the same length as the syntenic reference segment and is less than 2bp
        STABLESOURCE = ref_dict[ANCHORSEGMENTID]["StableSource"]

        if SEGMENTSENSE == "+":
            #print("positive")
            START = int(ref_dict[ANCHORSEGMENTID]["StableOffset"]) + int(ref_dict[ANCHORSEGMENTID]["SegmentLength"]) + SEGMENTOFFSET
        elif SEGMENTSENSE == "-":
            #print("negative")
            START = int(ref_dict[ANCHORSEGMENTID]["StableOffset"]) + int(ref_dict[ANCHORSEGMENTID]["SegmentLength"]) + query_dict[SEGMENTID]["SegmentLength"] - SEGMENTOFFSET
        else:
            print("Something went wrong")
            START = "NA"

        STOP = START + 1
        REFCHECK = False
        ANCHORCHECK = bool(ANCHORTEST)
        LENGTHCHECK = LENGTHTEST
        SYNTENICLENGTHCHECK = SYNTENICLENGHTTEST

    elif not REFTEST and not ANCHORTEST:
        # Segment is a query segment without an anchor
        STABLESOURCE = "NA"
        START = "NA"
        STOP = "NA"
        REFCHECK = False
        ANCHORCHECK = bool(ANCHORTEST)
        LENGTHCHECK = LENGTHTEST
        SYNTENICLENGTHCHECK = SYNTENICLENGHTTEST

    elif not REFTEST and ANCHORTEST and not LENGTHTEST and SYNTENICLENGHTTEST:
        # Segment is a query segment with an anchor and is the same length as the syntenic reference segment and is greater than 2bp
        STABLESOURCE = ref_dict[ANCHORSEGMENTID]["StableSource"]

        if SEGMENTSENSE == "+":
            #print("positive")
            START = int(ref_dict[ANCHORSEGMENTID]["StableOffset"]) + int(ref_dict[ANCHORSEGMENTID]["SegmentLength"]) + SEGMENTOFFSET
        elif SEGMENTSENSE == "-":
            #print("negative")
            START = int(ref_dict[ANCHORSEGMENTID]["StableOffset"]) + int(ref_dict[ANCHORSEGMENTID]["SegmentLength"]) + query_dict[SEGMENTID]["SegmentLength"] - SEGMENTOFFSET
        else:
            print("Something went wrong")
            START = "NA"

        STOP = START + 1
        REFCHECK = False
        ANCHORCHECK = bool(ANCHORTEST)
        LENGTHCHECK = LENGTHTEST
        SYNTENICLENGTHCHECK = SYNTENICLENGHTTEST

    elif not REFTEST and ANCHORTEST and LENGTHTEST and not SYNTENICLENGHTTEST:
        # Segment is a query segment with an anchor and is not the same length as the syntenic reference segment and is less than 2bp
        STABLESOURCE = ref_dict[ANCHORSEGMENTID]["StableSource"]

        if SEGMENTSENSE == "+":
            #print("positive")
            START = int(ref_dict[ANCHORSEGMENTID]["StableOffset"]) + int(ref_dict[ANCHORSEGMENTID]["SegmentLength"]) + SEGMENTOFFSET
        elif SEGMENTSENSE == "-":
            #print("negative")
            START = int(ref_dict[ANCHORSEGMENTID]["StableOffset"]) + int(ref_dict[ANCHORSEGMENTID]["SegmentLength"]) + query_dict[SEGMENTID]["SegmentLength"] - SEGMENTOFFSET
        else:
            print("Something went wrong")
            START = "NA"

        STOP = START + 1
        REFCHECK = False
        ANCHORCHECK = bool(ANCHORTEST)
        LENGTHCHECK = LENGTHTEST
        SYNTENICLENGTHCHECK = SYNTENICLENGHTTEST

    elif not REFTEST and ANCHORTEST and not LENGTHTEST and not SYNTENICLENGHTTEST:
        # Segment is a query segment with an anchor and is not the same length as the syntenic reference segment and is greater than 2bp
        STABLESOURCE = ref_dict[ANCHORSEGMENTID]["StableSource"]

        if SEGMENTSENSE == "+":
            #print("positive")
            START = int(ref_dict[ANCHORSEGMENTID]["StableOffset"]) + int(ref_dict[ANCHORSEGMENTID]["SegmentLength"]) + SEGMENTOFFSET
        elif SEGMENTSENSE == "-":
            #print("negative")
            START = int(ref_dict[ANCHORSEGMENTID]["StableOffset"]) + int(ref_dict[ANCHORSEGMENTID]["SegmentLength"]) + query_dict[SEGMENTID]["SegmentLength"] - SEGMENTOFFSET
        else:
            print("Something went wrong")
            START = "NA"

        STOP = START + 1
        REFCHECK = False
        ANCHORCHECK = bool(ANCHORTEST)
        LENGTHCHECK = LENGTHTEST
        SYNTENICLENGTHCHECK = SYNTENICLENGHTTEST

    else:
        print("Something went wrong")
        FLAG = True

    # This section will compute conversionCode using the convertFlag function
    CONVERSIONCODE = convertFlagsToFlagCode(REFCHECK, ANCHORCHECK, LENGTHCHECK, SYNTENICLENGTHCHECK, DOUBLEANCHORTEST, QUERYCHECK, FLAG)
    return ([STABLESOURCE, START, STOP, CONTEXT, METHYLATEDFRACTION, SENSE, COVERAGE, CONVERSIONCODE, SEGMENTID, SEGMENTOFFSET])


### This section will open GraphMethyl file
# if the RefOnlyParam is set to True, then the line will be processed through pullRefOnlyCoords
# if the RefOnlyParam is set to False, then the line will be processed through pullAllCoords
# processed lines will be saved to the OutputFile in tab-delimited format

print("Starting conversion...")
with open(OutputFile, "w", buffering=1000000) as output:
    with open(GraphMethylFile, "r") as input:
        for line in input:
            if RefOnlyParam == "True":
                output.write("\t".join([str(i) for i in pullRefOnlyCoords(line)]) + "\n")
            else:
                output.write("\t".join([str(i) for i in pullAllCoords(line)]) + "\n")

output.close()
input.close()

# with open(GraphMethylFile) as f:
#     contents=f.read().splitlines()
#
# if RefOnlyParam == "True":
#     print("Starting conversion...")
#     print("Only return reference segments")
#     OUTPUT = np.vstack([np.vstack([pullRefOnlyCoords(ROW) for ROW in contents])])
#     print("Saving result to file...")
#     np.savetxt(OutputFile, OUTPUT, delimiter="\t", fmt="%s")
# else:
#     print("Starting conversion...")
#     print("All segments")
#     OUTPUT = np.vstack([np.vstack([pullAllCoords(ROW) for ROW in contents])])
#     print("Saving result to file...")
#     np.savetxt(OutputFile, OUTPUT, delimiter="\t", fmt="%s")



print("Complete!")