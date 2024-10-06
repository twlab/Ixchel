#!/bin/bash
#SBATCH --mem=80G
#SBATCH -n 1
#SBATCH -N 1


echo "Loading software..."

eval $( spack load --sh python@3 )
eval $( spack load --sh py-numpy/i7mcgz4 )
eval $( spack load --sh py-jsonpickle@1.4.1 )
eval $( spack load --sh py-dill@0.3.4 )

INPUTFILE=$1
SEGMENTSPICKLE=$2
QUERYSEGMENTPICK=$3
LINKSPICKLE=$4
OUTPUTFILE=$5
REFONLYPARAM=$6
UpstreamLinksArrayFile=$7
DownstreamLinksArrayFile=$8
DoubleAnchorsFile=$9

echo "Input file: $INPUTFILE"
echo "Segments pickle: $SEGMENTSPICKLE"
echo "Query segment pickle: $QUERYSEGMENTPICK"
echo "Links pickle: $LINKSPICKLE"
echo "Output file: $OUTPUTFILE"
echo "Ref only param: $REFONLYPARAM"
echo "Upstream links array file: "$UpstreamLinksArrayFile
echo "Downstream links array file: "$DownstreamLinksArrayFile
echo "Double anchors file: "$DoubleAnchorsFile

echo "Start python script..."
python3 /scratch/hllab/Juan/JuanMacias_General_Code/convertSegmentToRef.py $INPUTFILE $SEGMENTSPICKLE $QUERYSEGMENTPICK $LINKSPICKLE $OUTPUTFILE $REFONLYPARAM $UpstreamLinksArrayFile $DownstreamLinksArrayFile $DoubleAnchorsFile


echo "Python script complete!"
