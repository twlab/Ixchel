#!/bin/bash
#SBATCH --mem=100G
#SBATCH -n 1
#SBATCH -N 1


echo "Loading software..."

eval $( spack load --sh python@3 )
eval $( spack load --sh py-numpy/i7mcgz4 )
eval $( spack load --sh py-jsonpickle@1.4.1 )
eval $( spack load --sh py-dill@0.3.4 )

INPUTLINKSFILE=$1
OUTOUTPICKLE=$2
DOUBLEANCHORPICKLEFILE=$3

echo "Start conversion..."
python3 /scratch/hllab/Juan/General_Code/makeAnchorLinkHashPickle.py $INPUTLINKSFILE $OUTOUTPICKLE $DOUBLEANCHORPICKLEFILE

echo "Complete!"
