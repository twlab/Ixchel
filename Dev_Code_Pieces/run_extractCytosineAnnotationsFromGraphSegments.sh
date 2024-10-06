#!/bin/bash
#SBATCH --mem=100G
#SBATCH -n 1
#SBATCH -N 1


echo "Loading software..."
eval $( spack load --sh python@3 )
eval $( spack load --sh py-numpy/i7mcgz4 )

INPUTSEGMENTSFILE=$1
OUTOUTANNOTATIONSFILE=$2

echo "Start extraction..."
python3 /scratch/hllab/Juan/JuanMacias_General_Code/extractCytosineAnnotationsFromGraphSegments.py $INPUTSEGMENTSFILE $OUTOUTANNOTATIONSFILE

echo "Complete!"
