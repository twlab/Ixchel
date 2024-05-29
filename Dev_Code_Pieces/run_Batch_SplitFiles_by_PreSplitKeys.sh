#!/bin/bash
#SBATCH --mem=8G
#SBATCH -n 1
#SBATCH -N 1

read FILETOSPLIT KEYSFILESPREFIX KEYSDIRECTORY OUTPUTPREFIX< <( sed -n ${SLURM_ARRAY_TASK_ID}p $1 )

### This script is meant to take in files needing to be split.
# For every Keys file in the keys directory, it will it will take the subset of lines in file to spit which match the keys in the keys file
# and output them to a file using the output prefix and the suffix of the keys file name which is _ delimited

echo "File to split: "$FILETOSPLIT
echo "Keys files prefix: "$KEYSFILESPREFIX
echo "Keys directory: "$KEYSDIRECTORY
echo "Output prefix: "$OUTPUTPREFIX

echo "Starting split..."

# This section will loop through all the keys files in the keys directory
for keysfile in $KEYSDIRECTORY$KEYSFILESPREFIX*; do
    echo "Keys file: "$keysfile
    OUTPUTFILE=$OUTPUTPREFIX$(echo $keysfile | awk -F "_" '{print $NF}')
    echo "     Output file: "$OUTPUTFILE
    # This section will take the subset of lines in the file to split which match the keys in the keys file
    # and output them to a file using the output prefix and the suffix of the keys file name which is _ delimited
    grep -F -f $keysfile $FILETOSPLIT >$OUTPUTFILE
done

echo " Complete!"


