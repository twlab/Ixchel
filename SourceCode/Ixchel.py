#!/usr/bin/env python3

# import necessary modules
import argparse
import sys
import numpy as np
from collections import defaultdict
import pickle
import dill
import subprocess
import re
import os
import glob
import shutil
import sqlite3
from tqdm import tqdm  # Progress bar library

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
        for line in tqdm(f, desc="Extracting segments", unit=" lines"):
            if line.startswith('S'):
                f_out.write(line)
    print(f"Extraction complete. Saved to {output_file}")

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
        for line in tqdm(f, desc="Extracting annotations", unit=" lines"):
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

    # Use system commands to report some statistics
    segments_count = count_lines(input_file)
    annotations_count = count_lines(output_file)
    print(f"Number of segments in {input_file}: {segments_count}")
    print(f"Number of annotations in {output_file}: {annotations_count}")

# Function to split segments into reference and query only files
def split_segments(args):
    input_file = args.input
    ref_name = args.reference_name if args.reference_name else "GRCh38"
    ref_output = f"RefOnly.{input_file}"
    query_output = f"QueryOnly.{input_file}"
    ref_pattern = re.compile(f"SN:Z:{ref_name}")
    with open(input_file, 'r') as f, open(ref_output, 'w') as ref_out, open(query_output, 'w') as query_out:
        for line in tqdm(f, desc="Splitting segments", unit=" lines"):
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
    base = os.path.splitext(INPUTFILE)[0]
    OUTPUTFILE = f"{base}.pkl"

    print(f"... Input file: {INPUTFILE}")
    with open(INPUTFILE) as f:
        for line in tqdm(f, desc="Processing segments", unit=" lines"):
            L = line.strip().split()
            SEGMENTID = L[1]
            SEQUENCELENGTH = len(L[2])
            STABLESOURCE = L[3].split(":")[2]
            STABLEOFFSET = L[4].split(":")[2]

            if bed_dict[SEGMENTID]["StableSource"]:
                print("ERROR 1")
                print(bed_dict[SEGMENTID]["SegmentLength"])
            else:
                bed_dict[SEGMENTID]["StableSource"] = STABLESOURCE
                bed_dict[SEGMENTID]["SegmentLength"] = SEQUENCELENGTH
                bed_dict[SEGMENTID]["StableOffset"] = STABLEOFFSET

    print(f"... Saving to {OUTPUTFILE}")
    with open(OUTPUTFILE, "wb") as f:
        pickle.dump(bed_dict, f)
    print("Reference segment hash pickle complete")

def makeQuerySegmentHashPickle(args):
    bed_dict = rec_dd()

    print("Making query segment hash pickle")
    INPUTFILE = args.input
    base = os.path.splitext(INPUTFILE)[0]
    OUTPUTFILE = f"{base}.pkl"

    print(f"... Input file: {INPUTFILE}")
    with open(INPUTFILE) as f:
        for line in tqdm(f, desc="Processing segments", unit=" lines"):
            L = line.strip().split()
            SEGMENTID = L[1]
            SEQUENCELENGTH = len(L[2])
            STABLESOURCE = ""
            STABLEOFFSET = ""

            if bed_dict[SEGMENTID]["SegmentLength"]:
                print("ERROR 1")
                print(bed_dict[SEGMENTID]["SegmentLength"])
            else:
                bed_dict[SEGMENTID]["StableSource"] = STABLESOURCE
                bed_dict[SEGMENTID]["SegmentLength"] = SEQUENCELENGTH
                bed_dict[SEGMENTID]["StableOffset"] = STABLEOFFSET

    print(f"Saving to {OUTPUTFILE}")
    with open(OUTPUTFILE, "wb") as f:
        pickle.dump(bed_dict, f)
    print("Query segment hash pickle complete")

def extract_links(args):
    input_file = args.input
    output_file = f"Links.{input_file}"

    print(f"Extracting links from {input_file} to {output_file}")
    with open(input_file, 'r') as f, open(output_file, 'w') as f_out:
        for line in tqdm(f, desc="Extracting links", unit=" lines"):
            if line.startswith('L'):
                f_out.write(line)
    print(f"Links extraction complete. Saved to {output_file}")

def create_link_search_keys(refsegmentsfile):
    search_keys_file = "temp_link_search_keys.txt"
    seen = set()
    print(f"... Creating search keys from {refsegmentsfile}")
    if not os.path.exists(refsegmentsfile):
        print(f"Error: {refsegmentsfile} does not exist!")
        sys.exit(1)
    with open(refsegmentsfile, 'r') as infile:
        for line in tqdm(infile, desc="Generating search keys", unit=" lines"):
            if line.startswith('S'):
                parts = line.split('\t')
                if len(parts) > 1:
                    key = f"L\t{parts[1]}\t\n"
                    seen.add(key)
    print(f"... Writing search keys to {search_keys_file}")
    with open(search_keys_file, 'w') as outfile:
        for key in sorted(seen):
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
        keys = set(line.strip() + "\t" for line in keys_file)

    print(f"... Filtering links using search keys")
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in tqdm(infile, desc="Filtering links", unit=" lines"):
            search_term = "L\t" + line.split('\t')[1] + "\t"
            if search_term in keys:
                outfile.write(line)

    print(f"... Removing temporary search keys file")
    os.remove(search_keys_file)
    print(f"Links filtering complete. Saved to {output_file}")

def makeAnchorLinkHashPickle(args):
    print("Making anchor link hash pickle")
    bed_dict = rec_dd()
    doubleanchor_dict = rec_dd()

    INPUTFILE = args.input
    base = os.path.splitext(INPUTFILE)[0]
    OUTPUTFILE = base + ".pkl"
    DOUBLEANCHORFILE = "DoubleAnchored." + base + ".pkl"

    print(f"... Input file: {INPUTFILE}")
    with open(INPUTFILE) as f:
        for line in tqdm(f, desc="Processing anchor links", unit=" lines"):
            L = line.strip().split()
            REFSEGMENTID = L[1]
            QUERYSEGMENTID = L[3]
            REFSTRAND = L[2]
            QUERYSTRAND = L[4]
            TAG = L[5]

            if bed_dict[QUERYSEGMENTID]["RefSegmentID"]:
                doubleanchor_dict[REFSEGMENTID]["QuerySegmentID"] = QUERYSEGMENTID
                bed_dict[QUERYSEGMENTID]["RefSegmentID"] = REFSEGMENTID
                bed_dict[QUERYSEGMENTID]["RefStrand"] = REFSTRAND
                bed_dict[QUERYSEGMENTID]["QueryStrand"] = QUERYSTRAND
                bed_dict[QUERYSEGMENTID]["Tag"] = TAG
            else:
                bed_dict[QUERYSEGMENTID]["RefSegmentID"] = REFSEGMENTID
                bed_dict[QUERYSEGMENTID]["RefStrand"] = REFSTRAND
                bed_dict[QUERYSEGMENTID]["QueryStrand"] = QUERYSTRAND
                bed_dict[QUERYSEGMENTID]["Tag"] = TAG

    print(f"... Saving to {OUTPUTFILE}")
    with open(OUTPUTFILE, "wb") as f:
        pickle.dump(bed_dict, f)
    print(f"... Saving double anchor file to {DOUBLEANCHORFILE}")
    with open(DOUBLEANCHORFILE, "wb") as f:
        pickle.dump(doubleanchor_dict, f)
    print("Anchor link hash pickle complete")

def makeLinkArrayPickles(args):
    ReferenceSegmentsPickle = args.ReferenceSegmentsPickle
    LinksPickle = args.FilteredLinksPickle
    UpstreamOutputFile = "UpstreamArray." + ReferenceSegmentsPickle
    DownstreamOutputFile = "DownstreamArray." + ReferenceSegmentsPickle
    print(f"... Input files: {ReferenceSegmentsPickle} and {LinksPickle}")
    print(f"... Output files: {UpstreamOutputFile} and {DownstreamOutputFile}")

    ref_dict = {}
    link_dict = {}

    print("Adding links segments...")
    with open(LinksPickle, 'rb') as f:
        link_dict = pickle.load(f)

    print("Adding reference segments...")
    with open(ReferenceSegmentsPickle, 'rb') as f:
        ref_dict = pickle.load(f)

    print("Building arrays of upstream and downstream keys from link_dict...")
    upstreamkeyarray = []
    downstreamkeyarray = []

    for key, value in tqdm(link_dict.items(), desc="Building link arrays", unit=" segments"):
        if ref_dict[value["RefSegmentID"]] and ref_dict[key]:
            downstreamkeyarray.append(key)
            upstreamkeyarray.append(value["RefSegmentID"])

    print(f"Saving to {UpstreamOutputFile}")
    with open(UpstreamOutputFile, "wb") as f:
        pickle.dump(upstreamkeyarray, f)

    print(f"Saving to {DownstreamOutputFile}")
    with open(DownstreamOutputFile, "wb") as f:
        pickle.dump(downstreamkeyarray, f)
    print("Link array pickles complete")

def split_annotations_file(args):
    print("Splitting annotations file")
    input_file = args.input
    if not os.path.exists("split_annotations"):
        os.mkdir("split_annotations")
    base = os.path.splitext(input_file)[0]
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

            for line in tqdm(file, desc="Splitting annotations", unit=" lines"):
                if count % lines_per_chunk == 0:
                    if current_file:
                        current_file.close()
                    current_file = open(f"{output_prefix}_{file_number:05d}", 'w')
                    file_number += 1
                current_file.write(line)
                count += 1

            if current_file:
                current_file.close()

        print("... Split complete!")

    except Exception as e:
        print(f"An error occurred!: {e}")

def precompute_conversion(args):
    print("Precomputing conversion to " + args.AnnotationFile + ".converted")

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

    print("... Adding query segments")
    with open(QuerySegmentsPickle, 'rb') as f:
        query_dict = pickle.load(f)

    print("... Adding links segments")
    with open(LinksPickle, 'rb') as f:
        link_dict = pickle.load(f)

    print("... Adding double anchor segments")
    with open(DoubleAnchorFile, 'rb') as f:
        doubleanchor_dict = pickle.load(f)

    def pullRefOnlyCoords(line):
        L = line.strip().split()
        SEGMENTID = L[0]
        SEGMENTOFFSET = int(L[1])
        SENSE = L[2]
        CONTEXT = L[3]
        COVERAGE = int(L[6])
        METHYLATEDFRACTION = L[7]
        if ref_dict[SEGMENTID]:
            STABLESOURCE = ref_dict[SEGMENTID]["StableSource"]
            START = int(ref_dict[SEGMENTID]["StableOffset"]) + SEGMENTOFFSET
            STOP = START + 1
            REFCHECK = True
            return ([STABLESOURCE, START, STOP, CONTEXT, METHYLATEDFRACTION, SENSE, COVERAGE, REFCHECK])
        else:
            return ([False, False, False, False, False, False, False, False])

    def convertFlagsToFlagCode(isReferenceFlag, hasAnchorFlag, lengthFlag, lengthMatchFlag, hasMultipleAnchorsFlag,
                               isQueryFlag, generalErrorFlag):
        flagCode = 0
        if isReferenceFlag:
            flagCode += 1
        if hasAnchorFlag:
            flagCode += 2
        if lengthFlag:
            flagCode += 4
        if lengthMatchFlag:
            flagCode += 8
        if hasMultipleAnchorsFlag:
            flagCode += 16
        if isQueryFlag:
            flagCode += 32
        if generalErrorFlag:
            flagCode += 64
        return flagCode

    print("... Adding upstream links array")
    with open(UpstreamOutputFile, 'rb') as f:
        upstreamkeyarray = pickle.load(f)

    print("... Adding downstream links array")
    with open(DownstreamOutputFile, 'rb') as f:
        downstreamkeyarray = pickle.load(f)

    def pullAllCoords(line):
        L = line.strip().split()
        SEGMENTID = L[0]
        SEGMENTOFFSET = int(L[1])
        SENSE = L[2]
        CONTEXT = L[3]
        COVERAGE = int(L[6])
        METHYLATEDFRACTION = L[7]

        REFTEST = ref_dict[SEGMENTID]
        QUERYTEST = query_dict[SEGMENTID]
        FLAG = False
        QUERYCHECK = bool(QUERYTEST)

        if not REFTEST and QUERYTEST:
            try:
                LENGTHTEST = query_dict[SEGMENTID]["SegmentLength"] <= 2
            except (TypeError, ValueError):
                LENGTHTEST = False
                FLAG = True

        elif REFTEST:
            LENGTHTEST = ref_dict[SEGMENTID]["SegmentLength"] <= 2
            DOUBLEANCHORTEST = False
        else:
            DOUBLEANCHORTEST = False
            LENGTHTEST = False
            FLAG = True

        ANCHORTEST = link_dict[SEGMENTID]

        if ANCHORTEST and not REFTEST and QUERYTEST:
            ANCHORSEGMENTID = link_dict[SEGMENTID]["RefSegmentID"]
            SEGMENTSENSE = link_dict[SEGMENTID]["QueryStrand"]
            if doubleanchor_dict[ANCHORSEGMENTID]["QuerySegmentID"]:
                SYNTENICLENGHTTEST = False
                FLAG = True
                DOUBLEANCHORTEST = True
            else:
                DOUBLEANCHORTEST = False
                try:
                    SyntenicReferenceSegmentID = downstreamkeyarray[upstreamkeyarray.index(ANCHORSEGMENTID)]
                    SYNTENICLENGHTTEST = ref_dict[SyntenicReferenceSegmentID]["SegmentLength"] == query_dict[SEGMENTID]["SegmentLength"]
                except ValueError:
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

        if REFTEST:
            STABLESOURCE = ref_dict[SEGMENTID]["StableSource"]
            START = int(ref_dict[SEGMENTID]["StableOffset"]) + SEGMENTOFFSET
            STOP = START + 1
            REFCHECK = True
            ANCHORCHECK = bool(ANCHORTEST)
            LENGTHCHECK = LENGTHTEST
            SYNTENICLENGTHCHECK = SYNTENICLENGHTTEST

        elif not REFTEST and ANCHORTEST and LENGTHTEST and SYNTENICLENGHTTEST:
            STABLESOURCE = ref_dict[ANCHORSEGMENTID]["StableSource"]

            if SEGMENTSENSE == "+":
                START = int(ref_dict[ANCHORSEGMENTID]["StableOffset"]) + int(
                    ref_dict[ANCHORSEGMENTID]["SegmentLength"]) + SEGMENTOFFSET
            elif SEGMENTSENSE == "-":
                START = int(ref_dict[ANCHORSEGMENTID]["StableOffset"]) + int(
                    ref_dict[ANCHORSEGMENTID]["SegmentLength"]) + query_dict[SEGMENTID]["SegmentLength"] - SEGMENTOFFSET
            else:
                START = "NA"

            STOP = START + 1
            REFCHECK = False
            ANCHORCHECK = bool(ANCHORTEST)
            LENGTHCHECK = LENGTHTEST
            SYNTENICLENGTHCHECK = SYNTENICLENGHTTEST

        elif not REFTEST and not ANCHORTEST:
            STABLESOURCE = "NA"
            START = "NA"
            STOP = "NA"
            REFCHECK = False
            ANCHORCHECK = bool(ANCHORTEST)
            LENGTHCHECK = LENGTHTEST
            SYNTENICLENGTHCHECK = SYNTENICLENGHTTEST

        elif not REFTEST and ANCHORTEST and not LENGTHTEST and SYNTENICLENGHTTEST:
            STABLESOURCE = ref_dict[ANCHORSEGMENTID]["StableSource"]

            if SEGMENTSENSE == "+":
                START = int(ref_dict[ANCHORSEGMENTID]["StableOffset"]) + int(
                    ref_dict[ANCHORSEGMENTID]["SegmentLength"]) + SEGMENTOFFSET
            elif SEGMENTSENSE == "-":
                START = int(ref_dict[ANCHORSEGMENTID]["StableOffset"]) + int(
                    ref_dict[ANCHORSEGMENTID]["SegmentLength"]) + query_dict[SEGMENTID]["SegmentLength"] - SEGMENTOFFSET
            else:
                START = "NA"

            STOP = START + 1
            REFCHECK = False
            ANCHORCHECK = bool(ANCHORTEST)
            LENGTHCHECK = LENGTHTEST
            SYNTENICLENGTHCHECK = SYNTENICLENGHTTEST

        elif not REFTEST and ANCHORTEST and LENGTHTEST and not SYNTENICLENGHTTEST:
            STABLESOURCE = ref_dict[ANCHORSEGMENTID]["StableSource"]

            if SEGMENTSENSE == "+":
                START = int(ref_dict[ANCHORSEGMENTID]["StableOffset"]) + int(
                    ref_dict[ANCHORSEGMENTID]["SegmentLength"]) + SEGMENTOFFSET
            elif SEGMENTSENSE == "-":
                START = int(ref_dict[ANCHORSEGMENTID]["StableOffset"]) + int(
                    ref_dict[ANCHORSEGMENTID]["SegmentLength"]) + query_dict[SEGMENTID]["SegmentLength"] - SEGMENTOFFSET
            else:
                START = "NA"

            STOP = START + 1
            REFCHECK = False
            ANCHORCHECK = bool(ANCHORTEST)
            LENGTHCHECK = LENGTHTEST
            SYNTENICLENGTHCHECK = SYNTENICLENGHTTEST

        elif not REFTEST and ANCHORTEST and not LENGTHTEST and not SYNTENICLENGHTTEST:
            STABLESOURCE = ref_dict[ANCHORSEGMENTID]["StableSource"]

            if SEGMENTSENSE == "+":
                START = int(ref_dict[ANCHORSEGMENTID]["StableOffset"]) + int(
                    ref_dict[ANCHORSEGMENTID]["SegmentLength"]) + SEGMENTOFFSET
            elif SEGMENTSENSE == "-":
                START = int(ref_dict[ANCHORSEGMENTID]["StableOffset"]) + int(
                    ref_dict[ANCHORSEGMENTID]["SegmentLength"]) + query_dict[SEGMENTID]["SegmentLength"] - SEGMENTOFFSET
            else:
                START = "NA"

            STOP = START + 1
            REFCHECK = False
            ANCHORCHECK = bool(ANCHORTEST)
            LENGTHCHECK = LENGTHTEST
            SYNTENICLENGTHCHECK = SYNTENICLENGHTTEST

        else:
            FLAG = True

        CONVERSIONCODE = convertFlagsToFlagCode(REFCHECK, ANCHORCHECK, LENGTHCHECK, SYNTENICLENGTHCHECK, DOUBLEANCHORTEST, QUERYCHECK, FLAG)
        return ([STABLESOURCE, START, STOP, CONTEXT, METHYLATEDFRACTION, SENSE, COVERAGE, CONVERSIONCODE, SEGMENTID, SEGMENTOFFSET])

    print("... Starting conversion")
    with open(OutputFile, "w", buffering=1000000) as output:
        with open(AnnotationFile, "r") as input:
            for line in tqdm(input, desc="Precomputing conversions", unit=" lines"):
                if RefOnlyParam == "True":
                    output.write("\t".join([str(i) for i in pullRefOnlyCoords(line)]) + "\n")
                else:
                    output.write("\t".join([str(i) for i in pullAllCoords(line)]) + "\n")

    print("Conversion complete")

def SerializePrecomputedPositionsHash(args):
    INPUTPRECOMPUTEDFILE = args.precomputedfile
    base = os.path.splitext(INPUTPRECOMPUTEDFILE)[0]
    OUTPUTDBFILE = base + ".db"

    print('Reading in PreComputedPositionsFile')
    conn = sqlite3.connect(OUTPUTDBFILE)
    cursor = conn.cursor()

    cursor.execute('''
        CREATE TABLE IF NOT EXISTS conversion (
            segment_id TEXT,
            segment_offset TEXT,
            stable_source TEXT,
            start INTEGER,
            stop INTEGER,
            conversion_code INTEGER,
            PRIMARY KEY (segment_id, segment_offset)
        )
    ''')

    with open(INPUTPRECOMPUTEDFILE, 'r') as PreComputedPositionsFile:
        data = []
        for line in tqdm(PreComputedPositionsFile, desc="Serializing positions", unit=" lines"):
            line = line.strip().split('\t')
            # Handle cases where the value is 'NA'
            start = int(line[1]) if line[1] != 'NA' else None
            stop = int(line[2]) if line[2] != 'NA' else None
            conversion_code = int(line[7]) if line[7] != 'NA' else None

            # Append the entry even if some values are None (NULL in SQL)
            data.append((line[8], line[9], line[0], start, stop, conversion_code))

            if len(data) % 1000000 == 0:
                cursor.executemany('INSERT OR REPLACE INTO conversion VALUES (?, ?, ?, ?, ?, ?)', data)
                conn.commit()
                print(f'Processed {len(data)} lines so far.')
                data = []

        if data:
            cursor.executemany('INSERT OR REPLACE INTO conversion VALUES (?, ?, ?, ?, ?, ?)', data)
            conn.commit()

    print("Serialization complete")

def postprepcleanup(args):
    print("Cleaning up intermediate files")
    gfafile = args.input
    base = os.path.splitext(gfafile)[0]

    precomputed_raw_files = glob.glob(f"split_annotations/Annotations.Segments.{base}__*.converted")
    for file in precomputed_raw_files:
        os.remove(file)

    if not os.path.exists("precomputed"):
        os.mkdir("precomputed")
    precomputed_files = glob.glob(f"split_annotations/Annotations.Segments.{base}__*.pkl")
    for file in precomputed_files:
        shutil.move(file, "precomputed")

    files_to_remove = [
        f"Annotations.Segments.{base}.gfa",
        f"DoubleAnchored.FilteredLinks.Links.{base}.pkl",
        f"DownstreamArray.RefOnly.Segments.{base}.pkl",
        f"FilteredLinks.Links.{base}.gfa",
        f"FilteredLinks.Links.{base}.pkl",
        f"Links.{base}.gfa",
        f"QueryOnly.Segments.{base}.gfa",
        f"QueryOnly.Segments.{base}.pkl",
        f"RefOnly.Segments.{base}.gfa",
        f"RefOnly.Segments.{base}.pkl",
        f"Segments.{base}.gfa",
        f"UpstreamArray.RefOnly.Segments.{base}.pkl",
    ]
    for file in tqdm(files_to_remove, desc="Removing files", unit=" files"):
        if os.path.exists(file):
            os.remove(file)
    print("Cleanup complete")

def convertGraphMethylToMethylC(args):
    precomputedfile = args.precomputedfile
    if not os.path.exists(precomputedfile):
        print(f"Error: {precomputedfile} does not exist.")
        sys.exit(1)
    conn = sqlite3.connect(precomputedfile)
    cursor = conn.cursor()

    graphmethylFile = args.input
    if not os.path.exists(graphmethylFile):
        print(f"Error: {graphmethylFile} does not exist.")
        sys.exit(1)
    outputfile = os.path.splitext(graphmethylFile)[0] + ".methylc"
    print(f"Converting GraphMethyl file: {graphmethylFile} to MethylC format: {outputfile}")

    # Determine if we're running in a Slurm batch job or an interactive session
    interactive = os.isatty(sys.stdout.fileno())

    with open(graphmethylFile, 'r') as f, open(outputfile, 'w') as f_out:
        for line in tqdm(f, desc="Converting GraphMethyl to MethylC", unit=" lines", disable=not interactive):
            L = line.strip().split()
            segmentID = L[0]
            segmentOffset = L[1]
            strand = L[2]
            context = L[3]
            coverage = L[6]
            methylatedFraction = L[7]
            cursor.execute('''
                SELECT stable_source, start, stop, conversion_code 
                FROM conversion 
                WHERE segment_id = ? AND segment_offset = ?
            ''', (segmentID, segmentOffset))
            result = cursor.fetchone()

            if result:
                stableSource, start, stop, conversionCode = result
                f_out.write(f"{stableSource}\t{start}\t{stop}\t{context}\t{methylatedFraction}\t{strand}\t{coverage}\t{segmentID}:{segmentOffset}\t{conversionCode}\n")
            else:
                f_out.write(f"NA\tNA\tNA\tNA\tNA\tNA\tNA\t{segmentID}:{segmentOffset}\tNA\n")

    conn.close()
    print(f"Conversion to MethylC format complete. Saved to {outputfile}")

def convertConversionCodeSingle(args):
    conversionCode = args.conversioncode
    isReferenceFlag = False
    hasAnchorFlag = False
    lengthFlag = False
    lengthMatchFlag = False
    hasMultipleAnchorsFlag = False
    isQueryFlag = False
    generalErrorFlag = False
    if conversionCode >= 64:
        generalErrorFlag = True
        conversionCode -= 64
    if conversionCode >= 32:
        isQueryFlag = True
        conversionCode -= 32
    if conversionCode >= 16:
        hasMultipleAnchorsFlag = True
        conversionCode -= 16
    if conversionCode >= 8:
        lengthMatchFlag = True
        conversionCode -= 8
    if conversionCode >= 4:
        lengthFlag = True
        conversionCode -= 4
    if conversionCode >= 2:
        hasAnchorFlag = True
        conversionCode -= 2
    if conversionCode >= 1:
        isReferenceFlag = True
        conversionCode -= 1
    print("isReferenceFlag\thasAnchorFlag\tlengthFlag\tlengthMatchFlag\thasMultipleAnchorsFlag\tisQueryFlag\tgeneralErrorFlag")
    print(f"{isReferenceFlag}\t{hasAnchorFlag}\t{lengthFlag}\t{lengthMatchFlag}\t{hasMultipleAnchorsFlag}\t{isQueryFlag}\t{generalErrorFlag}")

def convertConversionCodes(args):
    methylCfile = args.input
    if not os.path.exists(methylCfile):
        print(f"Error: {methylCfile} does not exist.")
        sys.exit(1)
    outputfile = os.path.splitext(methylCfile)[0] + ".conversioncodes"
    print(f"Converting conversion codes from {methylCfile} to {outputfile}")
    with open(methylCfile, 'r') as f, open(outputfile, 'w') as f_out:
        f_out.write("chromosome\tstart\tstop\tconversionCode\tisReferenceFlag\thasAnchorFlag\tlengthFlag\tlengthMatchFlag\thasMultipleAnchorsFlag\tisQueryFlag\tgeneralErrorFlag\n")
        for line in tqdm(f, desc="Converting conversion codes", unit=" lines"):
            L = line.strip().split()
            chromosome = L[0]
            start = L[1]
            stop = L[2]
            conversionCode = int(L[8])
            isReferenceFlag = False
            hasAnchorFlag = False
            lengthFlag = False
            lengthMatchFlag = False
            hasMultipleAnchorsFlag = False
            isQueryFlag = False
            generalErrorFlag = False
            if conversionCode >= 64:
                generalErrorFlag = True
                conversionCode -= 64
            if conversionCode >= 32:
                isQueryFlag = True
                conversionCode -= 32
            if conversionCode >= 16:
                hasMultipleAnchorsFlag = True
                conversionCode -= 16
            if conversionCode >= 8:
                lengthMatchFlag = True
                conversionCode -= 8
            if conversionCode >= 4:
                lengthFlag = True
                conversionCode -= 4
            if conversionCode >= 2:
                hasAnchorFlag = True
                conversionCode -= 2
            if conversionCode >= 1:
                isReferenceFlag = True
                conversionCode -= 1
            f_out.write(f"{chromosome}\t{start}\t{stop}\t{L[8]}\t{isReferenceFlag}\t{hasAnchorFlag}\t{lengthFlag}\t{lengthMatchFlag}\t{hasMultipleAnchorsFlag}\t{isQueryFlag}\t{generalErrorFlag}\n")
    print(f"Conversion codes processing complete. Saved to {outputfile}")

def prepareGraphFiles(args):
    gfafile = args.input

    ### Extract Segments
    extract_segments(args)

    ### Extract Annotations
    args.input = f"{gfafile}"
    extract_cytosine_annotations(args)

    ### Split Segments
    args.input = f"Segments.{gfafile}"
    print(f"\nSplitting segments from {args.input}...")
    split_segments(args)

    ### Serialize Query Segments
    args.input = f"QueryOnly.Segments.{gfafile}"
    print(f"\nSerializing Query Segments from {args.input}...")
    makeQuerySegmentHashPickle(args)

    ### Serialize Reference Segments
    args.input = f"RefOnly.Segments.{gfafile}"
    print(f"\nSerializing Reference Segments from {args.input}...")
    makeRefSegmentHashPickle(args)

    ### Extract Links
    args.input = f"{gfafile}"
    extract_links(args)

    ### Filter Links
    args.input = f"Links.{gfafile}"
    args.refsegmentsfile = f"RefOnly.Segments.{gfafile}"
    filter_links(args)

    ### Serialize Anchor Links
    args.input = f"FilteredLinks.Links.{gfafile}"
    print(f"\nSerializing Anchor Links from {args.input}...")
    makeAnchorLinkHashPickle(args)

    ### Make up and down stream links
    gfafileBase = os.path.splitext(gfafile)[0]
    args.ReferenceSegmentsPickle = f"RefOnly.Segments.{gfafileBase}.pkl"
    args.FilteredLinksPickle = f"FilteredLinks.Links.{gfafileBase}.pkl"
    print(f"\nMaking link array pickles from {args.ReferenceSegmentsPickle} and {args.FilteredLinksPickle}")
    makeLinkArrayPickles(args)

    ### Split Annotations
    args.input = f"Annotations.Segments.{gfafile}"
    split_annotations_file(args)

    ### Precompute conversions - all
    gfafileBase = os.path.splitext(gfafile)[0]
    args.AnnotationFile = f"Annotations.Segments.{gfafile}"
    args.ReferenceSegmentsPickle = f"RefOnly.Segments.{gfafileBase}.pkl"
    args.QuerySegmentsPickle = f"QueryOnly.Segments.{gfafileBase}.pkl"
    args.LinksPickle = f"FilteredLinks.Links.{gfafileBase}.pkl"
    args.RefOnlyParam = "False"
    args.UpstreamOutputFile = f"UpstreamArray.RefOnly.Segments.{gfafileBase}.pkl"
    args.DownstreamOutputFile = f"DownstreamArray.RefOnly.Segments.{gfafileBase}.pkl"
    args.DoubleAnchorFile = f"DoubleAnchored.FilteredLinks.Links.{gfafileBase}.pkl"
    print(f"\nPrecomputing conversion from {args.AnnotationFile} using {args.ReferenceSegmentsPickle}, {args.QuerySegmentsPickle}, {args.LinksPickle}")
    print(f"UpstreamOutputFile: {args.UpstreamOutputFile}")
    print(f"DownstreamOutputFile: {args.DownstreamOutputFile}")
    print(f"DoubleAnchorFile: {args.DoubleAnchorFile}")
    precompute_conversion(args)

    ### Build and serialize precomputed conversion dictionary
    args.precomputedfile = f"Annotations.Segments.{gfafile}.converted"
    print(f"\nSerializing precomputed conversion dictionary from {args.precomputedfile}")
    SerializePrecomputedPositionsHash(args)

    ### Post prep clean up function
    args.input = f"{gfafile}"
    print(f"\nCleaning up intermediate files from {args.input}")
    postprepcleanup(args)

def main():
    parser = argparse.ArgumentParser(description="Ixchel Tool for processing genome graphs")
    subparsers = parser.add_subparsers(dest='command', help='sub-command help')

    parser_segments = subparsers.add_parser('extract_segments', help='extract segment lines from a graph')
    parser_segments.add_argument('input', type=str, help='GFA file to extract segments from')
    parser_segments.set_defaults(func=extract_segments)

    parser_extract = subparsers.add_parser('extract_annotations', help='extract cytosine annotations from a graph')
    parser_extract.add_argument('input', type=str, help='Segments file to extract annotations from')
    parser_extract.set_defaults(func=extract_cytosine_annotations)

    parser_split = subparsers.add_parser('split_segments', help='split segments into RefOnly and QueryOnly files based on reference name')
    parser_split.add_argument('input', type=str, help='Segments file to split')
    parser_split.add_argument('--reference_name', type=str, help='reference name to filter by, default is GRCh38', default='GRCh38')
    parser_split.set_defaults(func=split_segments)

    parser_pickle = subparsers.add_parser('makeRefSegmentHashPickle', help='make a reference segment hash pickle')
    parser_pickle.add_argument('input', type=str, help='Segments file to serialize')
    parser_pickle.set_defaults(func=makeRefSegmentHashPickle)

    parser_pickle = subparsers.add_parser('makeQuerySegmentHashPickle', help='make a query segment hash pickle')
    parser_pickle.add_argument('input', type=str, help='Segments file to serialize')
    parser_pickle.set_defaults(func=makeQuerySegmentHashPickle)

    parser_links = subparsers.add_parser('extract_links', help='extract link lines from a graph')
    parser_links.add_argument('input', type=str, help='GFA file to extract links from')
    parser_links.set_defaults(func=extract_links)

    parser_links = subparsers.add_parser('filter_links', help='Filter links to set where the source is a reference segment')
    parser_links.add_argument('input', type=str, help='Links file to filter')
    parser_links.add_argument('refsegmentsfile', type=str, help='Reference segments file to use for filtering links')
    parser_links.set_defaults(func=filter_links)

    parser_pickle = subparsers.add_parser('makeAnchorLinkHashPickle', help='make an anchor link hash pickle')
    parser_pickle.add_argument('input', type=str, help='Filtered reference as source links file to serialize')
    parser_pickle.set_defaults(func=makeAnchorLinkHashPickle)

    parser_pickle = subparsers.add_parser('makeLinkArrayPickles', help='make link array pickles')
    parser_pickle.add_argument('ReferenceSegmentsPickle', type=str, help='Reference segments pickle file')
    parser_pickle.add_argument('FilteredLinksPickle', type=str, help='Filtered links pickle file')
    parser_pickle.set_defaults(func=makeLinkArrayPickles)

    parser_split = subparsers.add_parser('split_annotations_file', help='split annotations file into chunks')
    parser_split.add_argument('input', type=str, help='Annotations file to split')
    parser_split.add_argument('--lines_per_chunk', type=int, help='number of lines per chunk', default=100000)
    parser_split.set_defaults(func=split_annotations_file)

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

    parser_pickle = subparsers.add_parser('SerializePrecomputedPositionsHash', help='serialize precomputed positions hash')
    parser_pickle.add_argument('precomputedfile', type=str, help='Precomputed positions file to serialize')
    parser_pickle.set_defaults(func=SerializePrecomputedPositionsHash)

    parser_cleanup = subparsers.add_parser('post_prep_cleanup', help='post prep cleanup')
    parser_cleanup.add_argument('input', type=str, help='GFA file to clean up')
    parser_cleanup.set_defaults(func=postprepcleanup)

    parser_prepareGraphFiles = subparsers.add_parser('prepareGraphFiles', help='prepare graph files')
    parser_prepareGraphFiles.add_argument('input', type=str, help='GFA file to prepare')
    parser_prepareGraphFiles.add_argument('--reference_name', type=str, help='reference name to filter by, default is GRCh38', default='GRCh38')
    parser_prepareGraphFiles.add_argument('--lines_per_chunk', type=int, help='number of lines per chunk', default=100000)
    parser_prepareGraphFiles.set_defaults(func=prepareGraphFiles)

    parser_convertGraphMethylToMethylC = subparsers.add_parser('convertGraphMethylToMethylC', help='convert .graph.methyl file to .methylC')
    parser_convertGraphMethylToMethylC.add_argument('input', type=str, help='.graph.methyl file to convert')
    parser_convertGraphMethylToMethylC.add_argument('precomputedfile', type=str, help='serialized (.db) precomputed conversion dictionary file')
    parser_convertGraphMethylToMethylC.set_defaults(func=convertGraphMethylToMethylC)

    parser_convertConversionCode = subparsers.add_parser('convertConversionCodeSingle', help='convert single conversion code to flags\nthis is experimental still')
    parser_convertConversionCode.add_argument('conversioncode', type=int, help='conversion code to convert')
    parser_convertConversionCode.set_defaults(func=convertConversionCodeSingle)

    parser_convertConversionCodes = subparsers.add_parser('convertConversionCodes', help='convert all conversion codes in a .methylC file to flags\nthis is experimental and should not be used')
    parser_convertConversionCodes.add_argument('input', type=str, help='.methylC file to convert')
    parser_convertConversionCodes.set_defaults(func=convertConversionCodes)

    args = parser.parse_args()
    if not hasattr(args, 'func'):
        parser.print_help()
        sys.exit(1)

    args.func(args)

if __name__ == "__main__":
    main()
