#!/bin/bash
#SBATCH --mem=18G
#SBATCH -n 1
#SBATCH -N 1

read FILETOSPLIT KEYSFILE OUTPUTFILE < <( sed -n ${SLURM_ARRAY_TASK_ID}p $1 )

### This script is meant to take in files needing to be split.
# For the single File to split, it will take the subset of lines in file to spit which match the keys in the keys file
# and output them to a file using the output file name

echo "File to split: "$FILETOSPLIT
echo "Keys file to split by: "$KEYSFILE
echo "Output file: "$OUTPUTFILE

# this section will create a temporary file and store the first two columns of the input file in it (tab delimited)
# then it will pull the subset of lines in the keys file which match completely the keys in the temporary file
# these lines will be saved to a temporary keys file
# then the temporary keys file will be used to pull the subset of lines from the input file and save them to the output file
# the temporary files will be deleted at the end of the script
# the temporary files will be named the same as the input and keys files with a large random number appended to the end and the .tmp extension

# generate random number

RANDOMNUMBER=$RANDOM.$RANDOM.$RANDOM.$RANDOM

# Generate temporary file names
echo "Generating temporary file names..."
TEMPINPUTFILE=$FILETOSPLIT.$RANDOMNUMBER.tmp
TEMPKEYSFILE=$OUTPUTFILE.$RANDOMNUMBER.keys.tmp
echo "Temporary input file: "$TEMPINPUTFILE
echo "Temporary keys file: "$TEMPKEYSFILE

# create temporary file
echo "Creating temporary file..."
cut -d $'\t' -f 1,2 $FILETOSPLIT >$TEMPINPUTFILE

# pull subset of lines from keys file which match the keys in the temporary file requiring a full word match
echo "Pulling subset of lines from keys file..."
grep -F -f $TEMPINPUTFILE -w $KEYSFILE >$TEMPKEYSFILE



# pull subset of lines from input file which match the keys in the temporary keys file
echo "Starting split..."
grep -F -f $TEMPKEYSFILE -w $FILETOSPLIT >$OUTPUTFILE

# delete temporary files
echo "Deleting temporary files..."
rm $TEMPINPUTFILE
rm $TEMPKEYSFILE

echo " Complete!"