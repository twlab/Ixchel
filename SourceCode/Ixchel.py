#!/usr/bin/env python3

# import necessary modules
import argparse
import sys
import numpy as np
from collections import defaultdict
import pickle
import dill
import subprocess
import re  # Importing the regular expressions library for pattern matching
import os

# Nested dictionary function
def rec_dd():
    return defaultdict(rec_dd)

# Define function to generate entry for each cytosine position with "+" sense
def generate_entry_cytosine(segment_id, position):
    return [segment_id, position, "+", "C", 0, 0, 0, 0]


# Define function to generate entry for each cytosine position with "-" sense
def generate_entry_cytosine_reverse(segment_id, position):
    return [segment_id, position, "-", "C", 0, 0, 0, 0]


# Function to extract segment lines and save them to a new file
def extract_segments(args):
    input_file = args.input
    output_file = f"Segments.{input_file}"

    print(f"Extracting segments from {input_file} to {output_file}")
    with open(input_file, 'r') as f, open(output_file, 'w') as f_out:
        for line in f:
            if line.startswith('S'):
                f_out.write(line)


# Use subprocess to call wc -l for counting lines
def count_lines(filename):
    result = subprocess.run(['wc', '-l', filename], stdout=subprocess.PIPE, text=True)
    return int(result.stdout.split()[0])


# Main function to extract cytosine annotations
def extract_cytosine_annotations(args):
    input_file = f"Segments.{args.input}"
    output_file = f"Annotations.{input_file}"

    print(f"Extracting cytosine annotations from {input_file} to {output_file}")
    with open(input_file, 'r') as f, open(output_file, "w", buffering=1000000) as f_out:
        for line in f:
            if line.startswith('S'):  # Only process segment lines
                segment_id, sequence = line.strip().split()[1:3]
                cytosine_positions = [i for i in range(len(sequence)) if sequence[i] == "C"]
                guanine_positions = [i for i in range(len(sequence)) if sequence[i] == "G"]

                # Write cytosine entries
                for pos in cytosine_positions:
                    entry = generate_entry_cytosine(segment_id, pos)
                    f_out.write("\t".join(map(str, entry)) + "\n")

                # Write guanine entries (reverse complement positions)
                for pos in guanine_positions:
                    entry = generate_entry_cytosine_reverse(segment_id, pos)
                    f_out.write("\t".join(map(str, entry)) + "\n")

    # use system commands to report some statistics
    ## How many segments are there?
    segments_count = count_lines(input_file)
    print(f"    Number of segments in {input_file}: {segments_count}")
    ## How many Cytosine annotations need to be converted?
    annotations_count = count_lines(output_file)
    print(f"    Number of annotations in {output_file}: {annotations_count}")


# New function to split segments into reference and query only files
def split_segments(args):
    input_file = args.input
    ref_name = args.reference_name if args.reference_name else "GRCh38"
    ref_output = f"RefOnly.{input_file}"
    query_output = f"QueryOnly.{input_file}"

    ref_pattern = re.compile(rf"SN:Z:{ref_name}\.\w+")

    with open(input_file, 'r') as f, open(ref_output, 'w') as ref_out, open(query_output, 'w') as query_out:
        for line in f:
            if ref_pattern.search(line):
                ref_out.write(line)
            else:
                query_out.write(line)

    print(f"Reference lines written to {ref_output}: {count_lines(ref_output)}")
    print(f"Query lines written to {query_output}: {count_lines(query_output)}")

def makeRefSegmentHashPickle(args):
    bed_dict = rec_dd()

    print("Making reference segment hash pickle")
    INPUTFILE = args.input
    base = os.path.splitext(INPUTFILE)[0]  # Removes the current extension
    OUTPUTFILE = f"{base}.pkl"

    print("... Input file: ")
    print(INPUTFILE)

    with open(INPUTFILE) as f:
        for line in f:
            L = line.strip().split()
            SEGMENTID = L[1]
            SEQUENCELENGTH = len(L[2])
            STABLESOURCE = L[3].split(":")[2]
            STABLEOFFSET = L[4].split(":")[2]

            # It looks like for the purposes of merging pickles we have to make the outer most layer ALIGNMENT
            if bed_dict[SEGMENTID]["StableSource"]:
                # Not first position for read
                print("ERROR 1")
                print(bed_dict[SEGMENTID]["SegmentLength"])
            else:
                # First postion for read
                bed_dict[SEGMENTID]["StableSource"] = STABLESOURCE
                bed_dict[SEGMENTID]["SegmentLength"] = SEQUENCELENGTH
                bed_dict[SEGMENTID]["StableOffset"] = STABLEOFFSET

    # Save nested dictionary as pickle file
    print("... Saving to: ")
    print(OUTPUTFILE)

    f = open(OUTPUTFILE, "wb")
    pickle.dump(bed_dict, f)
    f.close()

def makeQuerySegmentHashPickle(args):
    bed_dict = rec_dd()

    print("Making query segment hash pickle")
    INPUTFILE = args.input
    base = os.path.splitext(INPUTFILE)[0]  # Removes the current extension
    OUTPUTFILE = f"{base}.pkl"

    print("... Input file: ")
    print(INPUTFILE)

    with open(INPUTFILE) as f:
        for line in f:
            L = line.strip().split()
            SEGMENTID = L[1]
            SEQUENCELENGTH = len(L[2])
            STABLESOURCE = ""
            STABLEOFFSET = ""
            # STABLESOURCE = L[3].split(":")[2]
            # STABLEOFFSET = L[4].split(":")[2]

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

def extract_links(args):
    input_file = args.input
    output_file = f"Links.{input_file}"

    print(f"Extracting links from {input_file} to {output_file}")
    with open(input_file, 'r') as f, open(output_file, 'w') as f_out:
        for line in f:
            if line.startswith('L'):
                f_out.write(line)

def create_link_search_keys(refsegmentsfile):
    search_keys_file = "temp_link_search_keys.txt"
    seen = set()  # This set will automatically handle unique entries
    print(f"... Creating search keys from {refsegmentsfile}")
    with open(refsegmentsfile, 'r') as infile:
        for line in infile:
            if line.startswith('S'):  # To mimic 'S\t"$2"\t' -> We take lines starting with 'S', THIS IS INEFFICIENT!!! REMOVE THIS
                parts = line.split('\t')
                if len(parts) > 1:
                    key = f"L\t{parts[1]}\t\n"  # Construct the key as per the awk command
                    seen.add(key)  # Add to set, which keeps entries unique
    print(f"... Writing search keys to {search_keys_file}")
    with open(search_keys_file, 'w') as outfile:
        for key in sorted(seen):  # Sort the set before writing
            outfile.write(key)

def filter_links(args):
    input_file = args.input
    refsegmentsfile = args.refsegmentsfile
    print(f"Filtering links from {input_file} using {refsegmentsfile}")
    output_file = f"FilteredLinks.{input_file}"
    print(f"... Generating search keys")
    create_link_search_keys(refsegmentsfile)
    search_keys_file = "temp_link_search_keys.txt"
    with open(search_keys_file, 'r') as keys_file:
        keys = set(line.strip() + "\t" for line in keys_file)  # Read all keys into a set

    print(f"... Filtering links using search keys")
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Construct the search term in the same format as keys are stored
            search_term = "L\t" + line.split('\t')[1] + "\t"
            if search_term in keys:
                outfile.write(line)  # Write to output if the search term is found in keys

    print(f"... Removing temporary search keys file")
    os.remove(search_keys_file)

def makeAnchorLinkHashPickle(args):
    print("Making anchor link hash pickle")
    bed_dict = rec_dd()
    doubleanchor_dict = rec_dd()

    INPUTFILE = args.input
    base = os.path.splitext(INPUTFILE)[0]  # Removes the current extension
    OUTPUTFILE = base + ".pkl"
    DOUBLEANCHORFILE = "DoubleAnchored." + base + ".pkl"

    print("... Input file:" + INPUTFILE)
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
                # print(["Double Anchor detected", QUERYSEGMENTID])
                doubleanchor_dict[REFSEGMENTID]["QuerySegmentID"] = QUERYSEGMENTID
                bed_dict[QUERYSEGMENTID]["RefSegmentID"] = REFSEGMENTID
                bed_dict[QUERYSEGMENTID]["RefStrand"] = REFSTRAND
                bed_dict[QUERYSEGMENTID]["QueryStrand"] = QUERYSTRAND
                bed_dict[QUERYSEGMENTID]["Tag"] = TAG
            else:
                # First position for segment
                bed_dict[QUERYSEGMENTID]["RefSegmentID"] = REFSEGMENTID
                bed_dict[QUERYSEGMENTID]["RefStrand"] = REFSTRAND
                bed_dict[QUERYSEGMENTID]["QueryStrand"] = QUERYSTRAND
                bed_dict[QUERYSEGMENTID]["Tag"] = TAG

    # Save nested dictionary as pickle file
    print("... Saving to:" + OUTPUTFILE)

    f = open(OUTPUTFILE, "wb")
    pickle.dump(bed_dict, f)
    f.close()

    print("... Saving double anchor file to:" + DOUBLEANCHORFILE)

    f = open(DOUBLEANCHORFILE, "wb")
    pickle.dump(doubleanchor_dict, f)
    f.close()

def makeLinkArrayPickles(args):
    print("Making link array pickles")
    ReferenceSegmentsPickle = args.ReferenceSegmentsPickle
    LinksPickle = args.FilteredLinksPickle
    UpstreamOutputFile = "UpstreamArray." + ReferenceSegmentsPickle
    DownstreamOutputFile = "DownstreamArray." + ReferenceSegmentsPickle
    print("... Input files:")
    print(ReferenceSegmentsPickle)
    print(LinksPickle)
    print("... Output files:")
    print(UpstreamOutputFile)
    print(DownstreamOutputFile)

    ref_dict = {}
    link_dict = {}

    print("Adding links segments...")
    with open(LinksPickle, 'rb') as f:
        link_dict = pickle.load(f)
    f.close()

    print("Adding reference segments...")
    with open(ReferenceSegmentsPickle, 'rb') as f:
        ref_dict = pickle.load(f)
    f.close()

    # Build arrays of upstream and downstream keys from link_dict
    print("Building arrays of upstream and downstream keys from link_dict...")
    upstreamkeyarray = []
    downstreamkeyarray = []

    for key, value in link_dict.items():
        # Filter out links that aren't reference segments to reference segments
        # print(ref_dict[value["RefSegmentID"]])
        # print(value["RefSegmentID"])
        # print(type(value["RefSegmentID"]))
        # print(key)
        # print(type(key))

        if ref_dict[value["RefSegmentID"]] and ref_dict[key]:
            downstreamkeyarray.append(key)
            upstreamkeyarray.append(value["RefSegmentID"])

        # 31611814
        # 31611815
        # if value["RefSegmentID"] == "31611815" and key == "31611814":
        #     #ref_dict[value["RefSegmentID"]] and ref_dict[key]:
        #     print("It matches")
        #     print(ref_dict[value["RefSegmentID"]])
        #     print(value["RefSegmentID"])
        #     print(type(value["RefSegmentID"]))
        #     print(repr(value["RefSegmentID"]))
        #     print(key)
        #     print(type(key))
        #     print(ref_dict[key])
        # print(upstreamkeyarray.index('31611815'))

    # Save nested dictionary as pickle file
    print("Saving to: ")

    print(UpstreamOutputFile)
    f = open(UpstreamOutputFile, "wb")
    pickle.dump(upstreamkeyarray, f)
    f.close()

    print(DownstreamOutputFile)
    f = open(DownstreamOutputFile, "wb")
    pickle.dump(downstreamkeyarray, f)
    f.close()


def split_annotations_file(args):
    print("Splitting annotations file")
    input_file = args.input
    if not os.path.exists("split_annotations"):
        os.mkdir("split_annotations")
    base = os.path.splitext(input_file)[0]  # Removes the current extension
    output_prefix = "split_annotations/" + base + "_"
    lines_per_chunk = args.lines_per_chunk
    print(f"... Input file: {input_file}")
    print(f"... Output prefix: {output_prefix}")
    print(f"... Lines per chunk: {lines_per_chunk}")

    try:
        with open(input_file, 'r') as file:
            count = 0
            file_number = 1
            current_file = None

            for line in file:
                if count % lines_per_chunk == 0:
                    if current_file:
                        current_file.close()
                    current_file = open(f"{output_prefix}_{file_number:05d}", 'w')  # 5 digits padding
                    file_number += 1
                current_file.write(line)
                count += 1

            if current_file:
                current_file.close()

        print("... Split complete!")

    except Exception as e:
        print(f"An error occurred!: {e}")


def precompute_conversion(args):
    print("Precomputing conversion of: " + args.AnnotationFile + " to " + args.AnnotationFile + ".converted")

    ref_dict = {}
    query_dict = {}
    link_dict = {}

    AnnotationFile = args.AnnotationFile
    ReferenceSegmentsPickle = args.ReferenceSegmentsPickle
    QuerySegmentsPickle = args.QuerySegmentsPickle
    LinksPickle = args.LinksPickle
    OutputFile = AnnotationFile + ".converted"
    RefOnlyParam = args.RefOnlyParam
    UpstreamOutputFile = args.UpstreamOutputFile
    DownstreamOutputFile = args.DownstreamOutputFile
    DoubleAnchorFile = args.DoubleAnchorFile

    print("... Adding reference segments")
    with open(ReferenceSegmentsPickle, 'rb') as f:
        ref_dict = pickle.load(f)
    f.close()

    print("... Adding query segments")
    with open(QuerySegmentsPickle, 'rb') as f:
        query_dict = pickle.load(f)
    f.close()

    print("... Adding links segments")
    with open(LinksPickle, 'rb') as f:
        link_dict = pickle.load(f)
    f.close()

    print("... Adding double anchor segments")
    with open(DoubleAnchorFile, 'rb') as f:
        doubleanchor_dict = pickle.load(f)
    f.close()

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
            return ([STABLESOURCE, START, STOP, CONTEXT, METHYLATEDFRACTION, SENSE, COVERAGE, REFCHECK])
        else:
            print("Query Segment")
            return ([False, False, False, False, False, False, False, False])

    ### This section will define the function that will convert the flags into a flag code.
    def convertFlagsToFlagCode(isReferenceFlag, hasAnchorFlag, lengthFlag, lengthMatchFlag, hasMultipleAnchorsFlag,
                               isQueryFlag, generalErrorFlag):
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

    print("... Adding upstream links array")
    with open(UpstreamOutputFile, 'rb') as f:
        upstreamkeyarray = pickle.load(f)
    f.close()

    print("... Adding downstream links array")
    with open(DownstreamOutputFile, 'rb') as f:
        downstreamkeyarray = pickle.load(f)
    f.close()

    def pullAllCoords(line):
        L = line.strip().split()
        SEGMENTID = L[0]
        SEGMENTOFFSET = int(L[1])
        SENSE = L[2]
        CONTEXT = L[3]
        # UNMETHYLATED = int(L[4])
        # METHYLATED = int(L[5])
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
                # This tests whether the anchor segment is a double anchor. For some reason the test isn't working.
                SYNTENICLENGHTTEST = False
                FLAG = True
                DOUBLEANCHORTEST = True
            else:
                DOUBLEANCHORTEST = False
                # This should pull the reference segment that is downstream of the anchor segment
                try:
                    SyntenicReferenceSegmentID = downstreamkeyarray[upstreamkeyarray.index(ANCHORSEGMENTID)]
                    SYNTENICLENGHTTEST = ref_dict[SyntenicReferenceSegmentID]["SegmentLength"] == query_dict[SEGMENTID][
                        "SegmentLength"]
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
                # print("positive")
                START = int(ref_dict[ANCHORSEGMENTID]["StableOffset"]) + int(
                    ref_dict[ANCHORSEGMENTID]["SegmentLength"]) + SEGMENTOFFSET
            elif SEGMENTSENSE == "-":
                # print("negative")
                START = int(ref_dict[ANCHORSEGMENTID]["StableOffset"]) + int(
                    ref_dict[ANCHORSEGMENTID]["SegmentLength"]) + query_dict[SEGMENTID]["SegmentLength"] - SEGMENTOFFSET
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
                # print("positive")
                START = int(ref_dict[ANCHORSEGMENTID]["StableOffset"]) + int(
                    ref_dict[ANCHORSEGMENTID]["SegmentLength"]) + SEGMENTOFFSET
            elif SEGMENTSENSE == "-":
                # print("negative")
                START = int(ref_dict[ANCHORSEGMENTID]["StableOffset"]) + int(
                    ref_dict[ANCHORSEGMENTID]["SegmentLength"]) + query_dict[SEGMENTID]["SegmentLength"] - SEGMENTOFFSET
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
                # print("positive")
                START = int(ref_dict[ANCHORSEGMENTID]["StableOffset"]) + int(
                    ref_dict[ANCHORSEGMENTID]["SegmentLength"]) + SEGMENTOFFSET
            elif SEGMENTSENSE == "-":
                # print("negative")
                START = int(ref_dict[ANCHORSEGMENTID]["StableOffset"]) + int(
                    ref_dict[ANCHORSEGMENTID]["SegmentLength"]) + query_dict[SEGMENTID]["SegmentLength"] - SEGMENTOFFSET
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
                # print("positive")
                START = int(ref_dict[ANCHORSEGMENTID]["StableOffset"]) + int(
                    ref_dict[ANCHORSEGMENTID]["SegmentLength"]) + SEGMENTOFFSET
            elif SEGMENTSENSE == "-":
                # print("negative")
                START = int(ref_dict[ANCHORSEGMENTID]["StableOffset"]) + int(
                    ref_dict[ANCHORSEGMENTID]["SegmentLength"]) + query_dict[SEGMENTID]["SegmentLength"] - SEGMENTOFFSET
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
        CONVERSIONCODE = convertFlagsToFlagCode(REFCHECK, ANCHORCHECK, LENGTHCHECK, SYNTENICLENGTHCHECK,
                                                DOUBLEANCHORTEST, QUERYCHECK, FLAG)
        return ([STABLESOURCE, START, STOP, CONTEXT, METHYLATEDFRACTION, SENSE, COVERAGE, CONVERSIONCODE, SEGMENTID,
                 SEGMENTOFFSET])

    ### This section will open GraphMethyl file
    # if the RefOnlyParam is set to True, then the line will be processed through pullRefOnlyCoords
    # if the RefOnlyParam is set to False, then the line will be processed through pullAllCoords
    # processed lines will be saved to the OutputFile in tab-delimited format

    print("... Starting conversion")
    with open(OutputFile, "w", buffering=1000000) as output:
        with open(AnnotationFile, "r") as input:
            for line in input:
                if RefOnlyParam == "True":
                    output.write("\t".join([str(i) for i in pullRefOnlyCoords(line)]) + "\n")
                else:
                    output.write("\t".join([str(i) for i in pullAllCoords(line)]) + "\n")

    output.close()
    input.close()
    print("... Complete!")

def SerializePrecomputedPositionsHash(args):
    print("Serializing PrecomputedPositionsHash")
    INPUTPRECOMPUTEDFILE = args.precomputedfile
    base = os.path.splitext(INPUTPRECOMPUTEDFILE)[0]  # Removes the current extension
    OUTPUTPICKLEFILE = base + ".pkl"
    ### This section will import the PreComputedPositionsFile and create a dictionary of the data in the PreComputedPositionsFile named ConversionDictionary. It will print a message when it has processed 1% increments of the file.
    ConversionDictionary = defaultdict(list)
    print('... Reading in PreComputedPositionsFile')
    with open(INPUTPRECOMPUTEDFILE, 'r') as PreComputedPositionsFile:
        for line in PreComputedPositionsFile:
            line = line.strip().split('\t')
            ConversionDictionary[(line[8], line[9])] = (line[0], line[1], line[2], line[7])
            if len(ConversionDictionary) % 1000000 == 0:
                print('Processed ' + str(len(ConversionDictionary)) + ' lines so far.')

    PreComputedPositionsFile.close()

    print('... Saving ConversionDictionary as a pickle file')
    ### This section will pickle ConversionDictionary and save it as a pickle file with the name in PickleFile.
    with open(OUTPUTPICKLEFILE, 'wb') as PickleFile:
        pickle.dump(ConversionDictionary, PickleFile, protocol=pickle.HIGHEST_PROTOCOL)
    PickleFile.close()


def main():
    parser = argparse.ArgumentParser(description="Ixchel Tool for processing genome graphs")
    subparsers = parser.add_subparsers(dest='command', help='sub-command help')

    # Parser for extracting segments
    parser_segments = subparsers.add_parser('extract_segments', help='extract segment lines from a graph')
    parser_segments.add_argument('input', type=str, help='GFA file to extract segments from')
    parser_segments.set_defaults(func=extract_segments)

    # Parser for extracting annotations
    parser_extract = subparsers.add_parser('extract_annotations', help='extract cytosine annotations from a graph')
    parser_extract.add_argument('input', type=str, help='Segments file to extract annotations from')
    parser_extract.set_defaults(func=extract_cytosine_annotations)

    # Parser for splitting segments into RefOnly and QueryOnly
    parser_split = subparsers.add_parser('split_segments', help='split segments into RefOnly and QueryOnly files based on reference name')
    parser_split.add_argument('input', type=str, help='Segments file to split')
    parser_split.add_argument('--reference_name', type=str, help='reference name to filter by, default is GRCh38', default='GRCh38')
    parser_split.set_defaults(func=split_segments)

    # Parser for making a reference segment hash pickle
    parser_pickle = subparsers.add_parser('makeRefSegmentHashPickle', help='make a reference segment hash pickle')
    parser_pickle.add_argument('input', type=str, help='Segments file to serialize')
    parser_pickle.set_defaults(func=makeRefSegmentHashPickle)

    # Parser for making a query segment hash pickle
    parser_pickle = subparsers.add_parser('makeQuerySegmentHashPickle', help='make a query segment hash pickle')
    parser_pickle.add_argument('input', type=str, help='Segments file to serialize')
    parser_pickle.set_defaults(func=makeQuerySegmentHashPickle)

    # Parser for extracting links
    parser_links = subparsers.add_parser('extract_links', help='extract link lines from a graph')
    parser_links.add_argument('input', type=str, help='GFA file to extract links from')
    parser_links.set_defaults(func=extract_links)

    # Parser for creating RefSourceLinks, make link search keys, then use them to filter all links
    parser_links = subparsers.add_parser('filter_links', help='Filter links to set where the source is a reference segment')
    parser_links.add_argument('input', type=str, help='Links file to filter')
    parser_links.add_argument('refsegmentsfile', type=str, help='Reference segments file to use for filtering links')
    parser_links.set_defaults(func=filter_links)

    # Parser for making an anchor link hash pickle
    parser_pickle = subparsers.add_parser('makeAnchorLinkHashPickle', help='make an anchor link hash pickle')
    parser_pickle.add_argument('input', type=str, help='Filtered reference as source links file to serialize')
    parser_pickle.set_defaults(func=makeAnchorLinkHashPickle)

    # Parser for making link array pickles
    parser_pickle = subparsers.add_parser('makeLinkArrayPickles', help='make link array pickles')
    parser_pickle.add_argument('ReferenceSegmentsPickle', type=str, help='Reference segments pickle file')
    parser_pickle.add_argument('FilteredLinksPickle', type=str, help='Filtered links pickle file')
    parser_pickle.set_defaults(func=makeLinkArrayPickles)

    # Parser for splitting annotations file. lines_per_chunk should have a default value
    parser_split = subparsers.add_parser('split_annotations_file', help='split annotations file into chunks')
    parser_split.add_argument('input', type=str, help='Annotations file to split')
    parser_split.add_argument('--lines_per_chunk', type=int, help='number of lines per chunk', default=100000)
    parser_split.set_defaults(func=split_annotations_file)

    # Parser for precomputing conversion
    parser_precompute = subparsers.add_parser('precompute_conversion', help='precompute conversion')
    parser_precompute.add_argument('AnnotationFile', type=str, help='Annotations file to convert')
    parser_precompute.add_argument('ReferenceSegmentsPickle', type=str, help='Reference segments pickle file')
    parser_precompute.add_argument('QuerySegmentsPickle', type=str, help='Query segments pickle file')
    parser_precompute.add_argument('LinksPickle', type=str, help='Links pickle file')
    parser_precompute.add_argument('--RefOnlyParam', type=str, help='RefOnlyParam', default='False')
    parser_precompute.add_argument('UpstreamOutputFile', type=str, help='UpstreamOutputFile')
    parser_precompute.add_argument('DownstreamOutputFile', type=str, help='DownstreamOutputFile')
    parser_precompute.add_argument('DoubleAnchorFile', type=str, help='DoubleAnchorFile')
    parser_precompute.set_defaults(func=precompute_conversion)

    # Parser for serializing precomputed positions hash
    parser_pickle = subparsers.add_parser('SerializePrecomputedPositionsHash', help='serialize precomputed positions hash')
    parser_pickle.add_argument('precomputedfile', type=str, help='Precomputed positions file to serialize')
    parser_pickle.set_defaults(func=SerializePrecomputedPositionsHash)

    args = parser.parse_args()
    if not hasattr(args, 'func'):
        parser.print_help()
        sys.exit(1)

    args.func(args)



if __name__ == "__main__":
    main()

