#!/bin/bash
#SBATCH --mem=12G
#SBATCH -n 1
#SBATCH -N 1

read GRAPHMETHYLFILE INPUTPRECOMPUTEDPICKLE OUTPUTMETHYLCFILE < <( sed -n ${SLURM_ARRAY_TASK_ID}p $1 )

echo "Loading software..."
eval $( spack load --sh python@3 )
eval $( spack load --sh py-jsonpickle@1.4.1 )
eval $( spack load --sh py-dill@0.3.4 )

#GRAPHMETHYLFILE=$1
#INPUTPRECOMPUTEDPICKLE=$2
#OUTPUTMETHYLCFILE=$3

echo "GraphMethyl file: $GRAPHMETHYLFILE"
echo "Input precomputed pickle file: $INPUTPRECOMPUTEDPICKLE"
echo "Output methylC file: $OUTPUTMETHYLCFILE"

echo "Start python script..."
python3 /scratch/hllab/Juan/JuanMacias_General_Code/convert_GraphMethyl_To_MethylC_fromPrecomputedPickle.py $GRAPHMETHYLFILE $INPUTPRECOMPUTEDPICKLE $OUTPUTMETHYLCFILE
echo "Python script complete!"

