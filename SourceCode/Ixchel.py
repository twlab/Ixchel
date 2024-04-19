#!/usr/bin/env python3

# import necessary modules
import argparse
import sys
import numpy as np
import jsonpickle
import dill


# Define extraction functions
## Define function to generate entry for each cytosine position with "+" sense
def generate_entry_cytosine(segment_id, position):
    entry = [segment_id, position, "+", "C", 0, 0, 0, 0]
    return entry

## Define function to generate entry for each cytosine position with "-" sense
def generate_entry_cytosine_reverse(segment_id, position):
    entry = [segment_id, position, "-", "C", 0, 0, 0, 0]
    return entry

## Define function to determine whether each segment has a relevant cytosine position
def process_segment_check_for_cytosine(line):
    L = line.strip().split()
    print(L)
    SEQUENCE = L[2]
    # determine if there are in "C" patten in the SEQUENCE string
    if "C" in SEQUENCE:
        return True
    else:
        return False

## Define function to determine whether each segment has a relevant guanine position
def process_segment_check_for_guanine(line):
    L = line.strip().split()
    SEQUENCE = L[2]

    # determine if there are in "G" patten in the SEQUENCE string
    if "G" in SEQUENCE:
        return True
    else:
        return False

## Define function to process each segment for cytosine positions
def process_segment_cytosine(line):
    L = line.strip().split()

    SEGMENTID = L[1]
    SEQUENCE = L[2]

    # find all instances of "C" pattern in the SEQUENCE string and return the index of the first character
    # in the pattern. POSITIONS are 0-based
    cytosine_positions = [i for i in range(len(SEQUENCE)) if SEQUENCE.startswith("C", i)]

    # Loopthrough all cytosine positions if any and generate an entry for each and place in ENTRIES
    if cytosine_positions:
        return np.vstack([generate_entry_cytosine(SEGMENTID, POSITION) for POSITION in cytosine_positions])

## Define function to process each segment for guanine positions
def process_segment_guanine(line):
    L = line.strip().split()

    SEGMENTID = L[1]
    SEQUENCE = L[2]

    # find all instances of "G" pattern in the SEQUENCE string and return the index of the first character
    # in the pattern. POSITIONS are 0-based
    guanine_positions = [i for i in range(len(SEQUENCE)) if SEQUENCE.startswith("G", i)]

    # Loopthrough all guanine positions if any and generate an entry for each and place in ENTRIES
    if guanine_positions:
        return np.vstack([generate_entry_cytosine_reverse(SEGMENTID, POSITION) for POSITION in guanine_positions])

## Define main function
def extract_cytosine_annotations(args):
    input_file = args.input
    output_file = args.output

    print(f"Extracting cytosine annotations from {input_file} to {output_file}")
    # Open input file and process each line
    with open(input_file) as f:
        # For each line ine the file, determine whether the corresponding segment has any cytosine positions and if so process that segment to generate entries. Then determine whether the corresponding segment has any guanine positions and if so process that segment to generate entries. Then generate a tab delimited file from all the entries.
        with open(output_file, "w", buffering=1000000) as f1:
            for line in f:
                if process_segment_check_for_cytosine(line):
                    entries = process_segment_cytosine(line)
                    if entries is not None:
                        for entry in entries:
                            f1.write("\t".join(map(str, entry)) + "\n")
                if process_segment_check_for_guanine(line):
                    entries = process_segment_guanine(line)
                    if entries is not None:
                        for entry in entries:
                            f1.write("\t".join(map(str, entry)) + "\n")
        f1.close()
    f.close()





def main():
    parser = argparse.ArgumentParser(description="Ixchel Tool for processing genome graphs")
    subparsers = parser.add_subparsers(dest='command', help='sub-command help')

    # create the parser for the "extract_annotations" command
    parser_extract = subparsers.add_parser('extract_annotations', help='extract cytosine annotations from a graph')
    parser_extract.add_argument('input', type=str, help='input GFA file')
    parser_extract.add_argument('output', type=str, help='output annotations file')
    parser_extract.set_defaults(func=extract_cytosine_annotations)

    args = parser.parse_args()
    if args.command is None:
        parser.print_help()
        sys.exit(1)

    # call the default function
    args.func(args)


if __name__ == "__main__":
    main()
