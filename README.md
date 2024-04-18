# Annotation surjection with Ixchel
Ixchel is a genome-graph based tool intended to aid in the conversion of annotations in graph-coordinates to linear coordinates.
It is early in development and currently only has functionality designed for surjecting annotations (.graph.methyl) data for CpG sites.

### Overview of Ixchel process
```mermaid
flowchart TB
    Genomegraph["Genome-graph (.GFA)"] --> IxchelPrep["1) Ixchel --prepareGraphFiles"]
    IxchelPrep --> GraphFiles["Pre-computed conversion files"]
    GraphFiles --> IxchelSurject["2) Ixchel --surject --mode CpG"]
    Annotations["Annotations (.graph.methyl)"] --> IxchelSurject
    IxchelSurject --> SurjectedAnnotations["Surjected annotations (.methylC)"]
    SurjectedAnnotations --> IxchelInterpretCodes["3) Ixchel --interpretCodes"]
    IxchelInterpretCodes --> ExpandedCodes["Expanded conversion codes"]
    
    classDef default fill:#fff,stroke:#333,stroke-width:4px,color:black;
```

## Prepare the graph files
Ixchel requires a genome graph in GFA format. It was designed to work with graphs produced by minigraph-cactus and is compatible with the human pangenomes available [here](https://github.com/human-pangenomics/hpp_pangenome_resources).

## Run surjection

## Extract and interpret codes

