#!/bin/bash
#SBATCH --mem=12G
#SBATCH -n 1
#SBATCH -N 1

read INPUTGFA REFERENCE OUTPUTRGFA < <( sed -n ${SLURM_ARRAY_TASK_ID}p $1 )

# Load the required module
eval $( spack load --sh /pudl6n3 )

echo "Job to convert GFA to rGFA for $INPUTGFA"
echo "With reference: $REFERENCE"

# Convert GFA to rGFA
python3 /scratch/hllab/Juan/JuanMacias_General_Code/Pangenomic/convert_GFA_to_rGFA.py \
$INPUTGFA \
-r $REFERENCE \
-o $OUTPUTRGFA
