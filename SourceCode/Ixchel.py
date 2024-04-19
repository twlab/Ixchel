#!/usr/bin/env python3

# import necessary modules
import argparse
import sys
import numpy as np
import jsonpickle
import dill


def extract_cytosine_annotations(args):
    input_file = args.input
    output_file = args.output

    print(f"Extracting cytosine annotations from {input_file} to {output_file}")

    # Example function body
    with open(input_file, 'r') as f, open(output_file, 'w') as f_out:
        for line in f:
            if 'C' in line:  # Simplified check
                f_out.write(line)  # Simplified processing logic


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
