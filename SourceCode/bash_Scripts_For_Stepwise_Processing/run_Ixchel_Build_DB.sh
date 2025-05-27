#!/bin/bash
#SBATCH --mem=8G
#SBATCH -n 1
#SBATCH -N 1

# Load the necessary modules
## python@3.7.3
eval $( spack load --sh /si3fu6h )
eval $( spack load --sh py-numpy/i7mcgz4 )
eval $( spack load --sh py-jsonpickle@1.4.1 )
## py-dill@0.3.4
eval $( spack load --sh py-dill/syrzqvm )
## py-tqdm
eval $( spack load --sh py-tqdm/vi7e33q )

# print input argument
echo "Input Precomputed conversion file: $1"
echo "Output db: $2"

# Run the script
echo "Running the script"
python3 /scratch/hllab/Juan/Ixchel/SourceCode/Ixchel.py build_db $1 $2
echo "Done"
