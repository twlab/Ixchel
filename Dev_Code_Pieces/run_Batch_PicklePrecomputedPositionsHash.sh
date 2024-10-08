#!/bin/bash
#SBATCH --mem=8G
#SBATCH -n 1
#SBATCH -N 1

read INPUTPRECOMPUTEDFILE OUTPUTPICKLEFILE < <( sed -n ${SLURM_ARRAY_TASK_ID}p $1 )

echo "Loading software..."
eval $( spack load --sh python@3 )
eval $( spack load --sh py-jsonpickle@1.4.1 )
eval $( spack load --sh py-dill@0.3.4 )

#INPUTPRECOMPUTEDFILE=$1
#OUTPUTPICKLEFILE=$2

echo "Input precomputed file: $INPUTPRECOMPUTEDFILE"
echo "Output pickle file: $OUTPUTPICKLEFILE"

echo "Start python script..."
python3 -u /scratch/hllab/Juan/JuanMacias_General_Code/makePrecomputedConversionPickle.py $INPUTPRECOMPUTEDFILE $OUTPUTPICKLEFILE

echo "Python script complete!"
