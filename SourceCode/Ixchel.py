#!/usr/bin/env python3

# import necessary modules
import argparse
import sys
import numpy as np
import jsonpickle
import dill
import subprocess
import re  # Importing the regular expressions library for pattern matching

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
    print(f"\n")
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


def main():
    parser = argparse.ArgumentParser(description="Ixchel Tool for processing genome graphs")
    subparsers = parser.add_subparsers(dest='command', help='sub-command help')

    # Parser for extracting segments
    parser_segments = subparsers.add_parser('extract_segments', help='extract segment lines from a graph')
    parser_segments.add_argument('input', type=str, help='input GFA file')
    parser_segments.set_defaults(func=extract_segments)

    # Parser for extracting annotations
    parser_extract = subparsers.add_parser('extract_annotations', help='extract cytosine annotations from a graph')
    parser_extract.add_argument('input', type=str, help='input GFA Segments file')
    parser_extract.set_defaults(func=extract_cytosine_annotations)

    # Parser for splitting segments into RefOnly and QueryOnly
    parser_split = subparsers.add_parser('split_segments', help='split segments into RefOnly and QueryOnly files based on reference name')
    parser_split.add_argument('input', type=str, help='input Segments file')
    parser_split.add_argument('--reference_name', type=str, help='reference name to filter by, default is GRCh38', default='GRCh38')
    parser_split.set_defaults(func=split_segments)

    args = parser.parse_args()
    if not hasattr(args, 'func'):
        parser.print_help()
        sys.exit(1)

    args.func(args)


if __name__ == "__main__":
    main()