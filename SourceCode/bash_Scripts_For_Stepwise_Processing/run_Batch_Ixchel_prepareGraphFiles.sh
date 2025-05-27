#!/bin/bash
#SBATCH --mem=64G
#SBATCH -n 1
#SBATCH -N 1

read GFAGRAPH REFERENCE < <( sed -n ${SLURM_ARRAY_TASK_ID}p $1 )
# GFAGRAPH is the .gfa format grah to be processed
# REFERENCE is the reference genome within the minigraph-cactus graph (i.e. GRCh38)

# Load the necessary modules
## python@3.7.3
eval $( spack load --sh /si3fu6h )
eval $( spack load --sh py-numpy/i7mcgz4 )
eval $( spack load --sh py-jsonpickle@1.4.1 )
eval $( spack load --sh py-dill/aldouqp )
eval $( spack load --sh py-tqdm/nfnho45 )

# print input argument
echo "Input Prepare Graph file for: $GFAGRAPH"
echo "Reference genome: $REFERENCE"

# Run the script
echo "Running the script"
python3 /scratch/hllab/Juan/Ixchel/SourceCode/Ixchel.py prepareGraphFiles $GFAGRAPH --reference_name $REFERENCE
echo "Done"
