#!/bin/bash
#SBATCH --mem=90G
#SBATCH -n 1
#SBATCH -N 1

read INPUTFILE SEGMENTSPICKLE QUERYSEGMENTSPICKLE LINKSPICKLE UpstreamLinksArrayFile DownstreamLinksArrayFile DoubleAnchorsFile < <( sed -n ${SLURM_ARRAY_TASK_ID}p $1 )

echo "Loading software..."

echo "Input file: "$INPUTFILE
echo "Ref segments pickle file: "$SEGMENTSPICKLE
echo "Query segments pickle file: "$QUERYSEGMENTSPICKLE
echo "Links pickle: "$LINKSPICKLE
echo "Upstream links array file: "$UpstreamLinksArrayFile
echo "Downstream links array file: "$DownstreamLinksArrayFile
echo "Double anchors file: "$DoubleAnchorsFile

eval $( spack load --sh python@3 )
eval $( spack load --sh py-numpy/i7mcgz4 )
eval $( spack load --sh py-jsonpickle@1.4.1 )
eval $( spack load --sh py-dill@0.3.4 )

echo "Start conversion..."
python3 /scratch/hllab/Juan/Ixchel/SourceCode/Ixchel.py precompute_conversion $INPUTFILE $SEGMENTSPICKLE $QUERYSEGMENTSPICKLE $LINKSPICKLE $UpstreamLinksArrayFile $DownstreamLinksArrayFile $DoubleAnchorsFile

echo "Complete!"
