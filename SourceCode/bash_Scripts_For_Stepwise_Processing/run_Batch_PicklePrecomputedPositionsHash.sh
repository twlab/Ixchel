#!/bin/bash
#SBATCH --mem=8G
#SBATCH -n 1
#SBATCH -N 1

read INPUTPRECOMPUTEDFILE < <( sed -n ${SLURM_ARRAY_TASK_ID}p $1 )

echo "... Loading software"
eval $( spack load --sh /pudl6n3 )
eval $( spack load --sh py-jsonpickle@1.4.1 )
eval $( spack load --sh py-dill@0.3.4 )

echo "... Input precomputed file: $INPUTPRECOMPUTEDFILE"

echo "... Start script..."
python3 -u /scratch/hllab/Juan/Ixchel/SourceCode/Ixchel.py SerializePrecomputedPositionsHash $INPUTPRECOMPUTEDFILE
echo "Complete!"
