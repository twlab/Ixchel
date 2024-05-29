#!/bin/bash
#SBATCH --mem=16G
#SBATCH -n 1
#SBATCH -N 1

### This script is meant to take in a parameter file that contains the following:
# 1. The input .graph.methyl file to be sorted
# 2. The output file name
### It will then sort the input file by segment and offset, and output the sorted file to the output file name

read UNSORTEDFILE SORTEDFILE< <( sed -n ${SLURM_ARRAY_TASK_ID}p $1 )

echo "Input file: "$UNSORTEDFILE
echo "Output file: "$SORTEDFILE

echo "Starting sort..."
sort -k1,1 -k2,2n $UNSORTEDFILE >$SORTEDFILE

echo " Complete!"