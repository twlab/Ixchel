#!/bin/bash
#SBATCH --mem=4G
#SBATCH -n 1
#SBATCH -N 1

read CONVERTEDPOSTIONSFILE SEGMENTANDOFFSETKEYSFILE < <( sed -n ${SLURM_ARRAY_TASK_ID}p $1 )

echo "Precomputed graph position conversion file in .methylC format: $CONVERTEDPOSTIONSFILE"
echo "File to save sorted and unique segment-ID's pairs: $SEGMENTANDOFFSETKEYSFILE"

echo "Sorting and removing duplicates from segment-ID's pairs..."
awk -F "\t" '{print $9"\t"$10"\t"}' $CONVERTEDPOSTIONSFILE | sort | uniq >$SEGMENTANDOFFSETKEYSFILE

echo "Done!"
