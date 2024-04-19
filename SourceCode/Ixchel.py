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
    base = os.path.splitext(input_file)[0]  # Removes the current extension
    output_prefix = base + "_"
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

    args = parser.parse_args()
    if not hasattr(args, 'func'):
        parser.print_help()
        sys.exit(1)

    args.func(args)



if __name__ == "__main__":
    main()

