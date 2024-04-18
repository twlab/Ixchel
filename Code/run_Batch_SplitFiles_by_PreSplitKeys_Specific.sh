#!/bin/bash
#SBATCH --mem=4G
#SBATCH -n 1
#SBATCH -N 1

read FILETOSPLIT KEYSFILE OUTPUTFILE < <( sed -n ${SLURM_ARRAY_TASK_ID}p $1 )

### This script is meant to take in files needing to be split.
# For the single File to split, it will take the subset of lines in file to spit which match the keys in the keys file
# and output them to a file using the output file name

echo "File to split: "$FILETOSPLIT
echo "Keys file to split by: "$KEYSFILE
echo "Output file: "$OUTPUTFILE

echo "Starting split..."
grep -F -f $KEYSFILE $FILETOSPLIT >$OUTPUTFILE

echo " Complete!"
