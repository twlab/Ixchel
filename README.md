# Annotation surjection with Ixchel
Ixchel is a genome-graph based tool intended to aid in the conversion of annotations in graph-coordinates to linear coordinates.
It is early in development and currently only has functionality designed for surjecting annotations (.graph.methyl) data for CpG sites.
Currently, it only works with CpG sites, but the goal is to expand it to other types of annotations.
Given the tremendous size of genomes, the steps will likely need to be run individually on a cluster.

## Overview of Ixchel process
```mermaid
flowchart TB
    Genomegraph["Genome-graph (.GFA)"] --> IxchelPrep["1) Ixchel prepareGraphFiles"]
    IxchelPrep --> GraphFiles["Pre-computed conversion files"]
    GraphFiles --> IxchelSurject["2) Ixchel convertGraphMethylToMethylC"]
    Annotations["Annotations (.graph.methyl)"] --> IxchelSurject
    IxchelSurject --> SurjectedAnnotations["Surjected annotations (.methylC)"]
    SurjectedAnnotations --> IxchelInterpretCodes["3) Ixchel convertConversionCodes"]
    IxchelInterpretCodes --> ExpandedCodes["Expanded conversion codes"]
    
    classDef default fill:#fff,stroke:#333,stroke-width:4px,color:black;
```

## Run Ixchel end-to-end
Ixchel requires a genome graph in GFA format. It was designed to work with graphs produced by minigraph-cactus and is compatible with the human pangenomes available [here](https://github.com/human-pangenomics/hpp_pangenome_resources).
It is run in two steps:
1. Prepare the graph files
2. Convert the graph.methyl annotations to methylC format
```bash
cd /scratch/hllab/Juan/Ixchel_Dev_Tests/VGPlayGround
srun --mem=8000 --cpus-per-task=1 -J interactive -p interactive --pty /bin/bash -l
eval $( spack load --sh python@3.7.3 )
eval $( spack load --sh py-numpy/i7mcgz4 )
eval $( spack load --sh py-jsonpickle@1.4.1 )
eval $( spack load --sh py-dill@0.3.4 )
```
### Prepare the graph files
```bash
python3 /scratch/hllab/Juan/Ixchel/SourceCode/Ixchel.py prepareGraphFiles TestGraph.gfa
```
### Convert GraphMethyl to MethylC
```bash
python3 /scratch/hllab/Juan/Ixchel/SourceCode/Ixchel.py convertGraphMethylToMethylC Example.CG.graph.methyl Annotations.Segments.TestGraph.gfa.pkl
```
### Extract and interpret codes
Ixchel uses a set of codes to represent the context of segments in the graph.
Context that affects how precisely the segment can be surjected to linear coordinates.
```bash
python3 /scratch/hllab/Juan/Ixchel/SourceCode/Ixchel.py convertConversionCodes Example.CG.graph.methylc
```

# More details
```mermaid
flowchart TB
  %%=== CLI Dispatch ===%%
  subgraph CLI["Ixchel.py CLI"]
    direction TB
    ES(extract_segments)
    EA(extract_annotations)
    SS(split_segments)
    MRP(makeRefSegmentHashlePickle)
    MQP(makeQuerySegmentHashPickle)
    XL(extract_links)
    FL(filter_links)
    MA(makeAnchorLinkHashPickle)
    ML(makeLinkArrayPickles)
    SA(split_annotations_file)
    PC(precompute_conversion)
    SPH(SerializePrecomputedPositionsHash)
    PCU(postprepcleanup)
    CGM(convertGraphMethylToMethylC)
    CCS(convertConversionCodeSingle)
    CCC(convertConversionCodes)
    PGP(prepareGraphFiles)
    BDB(build_db)
    CMO(convert_methyl_optimized)
  end

  %%=== PrepareGraphFiles sequence ===%%
  subgraph Prepare["prepareGraphFiles ➔"]
    direction LR
    ES2(extract_segments) --> EA2(extract_annotations) --> SS2(split_segments)
    SS2 --> MQP2(makeQuerySegmentHashlePickle) --> MRP2(makeRefSegmentHashlePickle)
    MRP2 --> XL2(extract_links) --> FL2(filter_links)
    FL2 --> MA2(makeAnchorLinkHashlePickle) --> ML2(makeLinkArrayPickles)
    ML2 --> SA2(split_annotations_file)
  end
  PGP --> Prepare

  %%=== Precompute Conversion internals ===%%
  subgraph PreComp["precompute_conversion ➔"]
    direction TB
    PRC[pullRefOnlyCoords]  
    PAC[pullAllCoords]
    CFC[convertFlagsToFlagCode]
    PRC & PAC & CFC --> PC
  end

  %%=== Original Serialize vs Optimized Build ===%%
  SPH -->|row-by-row| SQL0[(conversion table)]
  BDB -->|bulk INSERT| SQL1[(conversion table)]

  %%=== Convert GraphMethyl ===%%
  subgraph ConvOLD["convertGraphMethylToMethylC (old)"]
    direction TB
    CGM --> SQL2["SELECT per record"]
    SQL2 --> CGM
  end

  subgraph ConvOPT["convert_methyl_optimized (new)"]
    direction TB
    CMO --> GroupSegs[group by segment]
    GroupSegs --> SQL3["SELECT per segment"]
    SQL3 --> Dict[mapping: offset→tuple]
    Dict --> CMO
  end

  %%=== CLI connections ===%%
  CLI --> ES
  CLI --> EA
  CLI --> SS
  CLI --> MRP
  CLI --> MQP
  CLI --> XL
  CLI --> FL
  CLI --> MA
  CLI --> ML
  CLI --> SA
  CLI --> PC
  CLI --> SPH
  CLI --> PCU
  CLI --> CGM
  CLI --> CCS
  CLI --> CCC
  CLI --> PGP
  CLI --> BDB
  CLI --> CMO
```

*For further support contact: juanfmacias[at]wustl.edu*