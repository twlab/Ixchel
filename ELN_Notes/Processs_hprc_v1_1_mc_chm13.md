# Overview
Here I will describe the proecss of preparing the `hprc_v1_1_mc_chm13` graph for surjection of methylGrapher data.
Ixchel pre-processing needs to be done per chromosome. The per-chromosome graph files are available [here](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms)

This is using `version v1.60.0 "Annicco` of vg

## Setup working directory
```bash
cd /scratch/hllab/Juan/Ixchel_Dev_Tests/

mkdir -p Processs_hprc_v1_1_mc_chm13
```

## Download per-chromosome graph files
```bash
cd /scratch/hllab/Juan/Ixchel_Dev_Tests/Processs_hprc_v1_1_mc_chm13

srun --mem=16000 --cpus-per-task=1 -J interactive -p interactive --pty /bin/bash -l

# Using latest awscli from spack
eval $( spack load --sh awscli )

aws s3 sync s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms ./vg_files --no-sign-request --exclude "*" --include "*.vg"

exit
```
```console
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr10.vg to vg_files/chr10.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr10.d9.vg to vg_files/chr10.d9.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr11.d9.vg to vg_files/chr11.d9.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr11.vg to vg_files/chr11.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr1.vg to vg_files/chr1.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr12.vg to vg_files/chr12.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr12.d9.vg to vg_files/chr12.d9.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr1.d9.vg to vg_files/chr1.d9.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr13.vg to vg_files/chr13.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr13.d9.vg to vg_files/chr13.d9.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr14.vg to vg_files/chr14.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr14.d9.vg to vg_files/chr14.d9.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr15.d9.vg to vg_files/chr15.d9.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr15.vg to vg_files/chr15.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr16.vg to vg_files/chr16.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr17.vg to vg_files/chr17.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr17.d9.vg to vg_files/chr17.d9.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr16.d9.vg to vg_files/chr16.d9.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr18.d9.vg to vg_files/chr18.d9.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr18.vg to vg_files/chr18.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr19.vg to vg_files/chr19.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr19.d9.vg to vg_files/chr19.d9.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr20.vg to vg_files/chr20.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr21.d9.vg to vg_files/chr21.d9.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr20.d9.vg to vg_files/chr20.d9.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr21.vg to vg_files/chr21.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr22.vg to vg_files/chr22.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr22.d9.vg to vg_files/chr22.d9.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr2.vg to vg_files/chr2.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr2.d9.vg to vg_files/chr2.d9.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr3.vg to vg_files/chr3.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr3.d9.vg to vg_files/chr3.d9.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr4.d9.vg to vg_files/chr4.d9.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr4.vg to vg_files/chr4.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr5.d9.vg to vg_files/chr5.d9.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr5.vg to vg_files/chr5.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr6.vg to vg_files/chr6.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr6.d9.vg to vg_files/chr6.d9.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr7.vg to vg_files/chr7.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr7.d9.vg to vg_files/chr7.d9.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chrEBV.d9.vg to vg_files/chrEBV.d9.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chrEBV.vg to vg_files/chrEBV.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chrM.d9.vg to vg_files/chrM.d9.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chrM.vg to vg_files/chrM.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chrOther.d9.vg to vg_files/chrOther.d9.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chrOther.vg to vg_files/chrOther.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr8.d9.vg to vg_files/chr8.d9.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr8.vg to vg_files/chr8.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chrY.d9.vg to vg_files/chrY.d9.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chrY.vg to vg_files/chrY.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chrX.d9.vg to vg_files/chrX.d9.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr9.d9.vg to vg_files/chr9.d9.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr9.vg to vg_files/chr9.vg
download: s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chrX.vg to vg_files/chrX.vg
```

## Convert chunks to gfa
```bash
cd /scratch/hllab/Juan/Ixchel_Dev_Tests/Processs_hprc_v1_1_mc_chm13

ls -lh vg_files/
```
```console
total 119G
-rw-rw-r-- 1 juanfmacias hllab 5.2G Jul  6  2023 chr1.d9.vg
-rw-rw-r-- 1 juanfmacias hllab 4.3G Jul  6  2023 chr1.vg
-rw-rw-r-- 1 juanfmacias hllab 2.9G Jul  6  2023 chr10.d9.vg
-rw-rw-r-- 1 juanfmacias hllab 2.5G Jul  6  2023 chr10.vg
-rw-rw-r-- 1 juanfmacias hllab 3.0G Jul  6  2023 chr11.d9.vg
-rw-rw-r-- 1 juanfmacias hllab 2.4G Jul  6  2023 chr11.vg
-rw-rw-r-- 1 juanfmacias hllab 2.9G Jul  6  2023 chr12.d9.vg
-rw-rw-r-- 1 juanfmacias hllab 2.3G Jul  6  2023 chr12.vg
-rw-rw-r-- 1 juanfmacias hllab 2.2G Jul  6  2023 chr13.d9.vg
-rw-rw-r-- 1 juanfmacias hllab 1.8G Jul  6  2023 chr13.vg
-rw-rw-r-- 1 juanfmacias hllab 2.1G Jul  6  2023 chr14.d9.vg
-rw-rw-r-- 1 juanfmacias hllab 1.7G Jul  6  2023 chr14.vg
-rw-rw-r-- 1 juanfmacias hllab 2.0G Jul  6  2023 chr15.d9.vg
-rw-rw-r-- 1 juanfmacias hllab 1.5G Jul  6  2023 chr15.vg
-rw-rw-r-- 1 juanfmacias hllab 2.4G Jul  6  2023 chr16.d9.vg
-rw-rw-r-- 1 juanfmacias hllab 2.0G Jul  6  2023 chr16.vg
-rw-rw-r-- 1 juanfmacias hllab 1.9G Jul  6  2023 chr17.d9.vg
-rw-rw-r-- 1 juanfmacias hllab 1.5G Jul  6  2023 chr17.vg
-rw-rw-r-- 1 juanfmacias hllab 1.7G Jul  6  2023 chr18.d9.vg
-rw-rw-r-- 1 juanfmacias hllab 1.5G Jul  6  2023 chr18.vg
-rw-rw-r-- 1 juanfmacias hllab 1.7G Jul  6  2023 chr19.d9.vg
-rw-rw-r-- 1 juanfmacias hllab 1.3G Jul  6  2023 chr19.vg
-rw-rw-r-- 1 juanfmacias hllab 5.3G Jul  6  2023 chr2.d9.vg
-rw-rw-r-- 1 juanfmacias hllab 4.2G Jul  6  2023 chr2.vg
-rw-rw-r-- 1 juanfmacias hllab 1.6G Jul  6  2023 chr20.d9.vg
-rw-rw-r-- 1 juanfmacias hllab 1.2G Jul  6  2023 chr20.vg
-rw-rw-r-- 1 juanfmacias hllab 1.1G Jul  6  2023 chr21.d9.vg
-rw-rw-r-- 1 juanfmacias hllab 918M Jul  6  2023 chr21.vg
-rw-rw-r-- 1 juanfmacias hllab 1.1G Jul  6  2023 chr22.d9.vg
-rw-rw-r-- 1 juanfmacias hllab 893M Jul  6  2023 chr22.vg
-rw-rw-r-- 1 juanfmacias hllab 4.1G Jul  6  2023 chr3.d9.vg
-rw-rw-r-- 1 juanfmacias hllab 3.3G Jul  6  2023 chr3.vg
-rw-rw-r-- 1 juanfmacias hllab 4.5G Jul  6  2023 chr4.d9.vg
-rw-rw-r-- 1 juanfmacias hllab 3.8G Jul  6  2023 chr4.vg
-rw-rw-r-- 1 juanfmacias hllab 3.8G Jul  6  2023 chr5.d9.vg
-rw-rw-r-- 1 juanfmacias hllab 3.5G Jul  6  2023 chr5.vg
-rw-rw-r-- 1 juanfmacias hllab 3.9G Jul  6  2023 chr6.d9.vg
-rw-rw-r-- 1 juanfmacias hllab 3.0G Jul  6  2023 chr6.vg
-rw-rw-r-- 1 juanfmacias hllab 3.8G Jul  6  2023 chr7.d9.vg
-rw-rw-r-- 1 juanfmacias hllab 3.2G Jul  6  2023 chr7.vg
-rw-rw-r-- 1 juanfmacias hllab 3.5G Jul  6  2023 chr8.d9.vg
-rw-rw-r-- 1 juanfmacias hllab 2.7G Jul  6  2023 chr8.vg
-rw-rw-r-- 1 juanfmacias hllab 3.4G Jul  6  2023 chr9.d9.vg
-rw-rw-r-- 1 juanfmacias hllab 2.7G Jul  6  2023 chr9.vg
-rw-rw-r-- 1 juanfmacias hllab 3.5M Jul  6  2023 chrEBV.d9.vg
-rw-rw-r-- 1 juanfmacias hllab 3.2M Jul  6  2023 chrEBV.vg
-rw-rw-r-- 1 juanfmacias hllab 518K Jul  6  2023 chrM.d9.vg
-rw-rw-r-- 1 juanfmacias hllab 453K Jul  6  2023 chrM.vg
-rw-rw-r-- 1 juanfmacias hllab  20M Jul  6  2023 chrOther.d9.vg
-rw-rw-r-- 1 juanfmacias hllab  17M Jul  6  2023 chrOther.vg
-rw-rw-r-- 1 juanfmacias hllab 2.1G Jul  6  2023 chrX.d9.vg
-rw-rw-r-- 1 juanfmacias hllab 1.8G Jul  6  2023 chrX.vg
-rw-rw-r-- 1 juanfmacias hllab 120M Jul  6  2023 chrY.d9.vg
-rw-rw-r-- 1 juanfmacias hllab 107M Jul  6  2023 chrY.vg
```
```bash
cd /scratch/hllab/Juan/Ixchel_Dev_Tests/Processs_hprc_v1_1_mc_chm13
srun --mem=32000 --cpus-per-task=20 -J interactive -p interactive --pty /bin/bash -l
/scratch/hllab/Juan/Latest_vg/vg convert -f vg_files/chr1.vg >chr1.gfa
/scratch/hllab/Juan/Latest_vg/vg convert -f vg_files/chr2.vg >chr2.gfa
/scratch/hllab/Juan/Latest_vg/vg convert -f vg_files/chr3.vg >chr3.gfa
/scratch/hllab/Juan/Latest_vg/vg convert -f vg_files/chr4.vg >chr4.gfa
/scratch/hllab/Juan/Latest_vg/vg convert -f vg_files/chr5.vg >chr5.gfa
/scratch/hllab/Juan/Latest_vg/vg convert -f vg_files/chr6.vg >chr6.gfa
/scratch/hllab/Juan/Latest_vg/vg convert -f vg_files/chr7.vg >chr7.gfa
/scratch/hllab/Juan/Latest_vg/vg convert -f vg_files/chr8.vg >chr8.gfa
/scratch/hllab/Juan/Latest_vg/vg convert -f vg_files/chr9.vg >chr9.gfa
/scratch/hllab/Juan/Latest_vg/vg convert -f vg_files/chr10.vg >chr10.gfa
/scratch/hllab/Juan/Latest_vg/vg convert -f vg_files/chr11.vg >chr11.gfa
/scratch/hllab/Juan/Latest_vg/vg convert -f vg_files/chr12.vg >chr12.gfa
/scratch/hllab/Juan/Latest_vg/vg convert -f vg_files/chr13.vg >chr13.gfa
/scratch/hllab/Juan/Latest_vg/vg convert -f vg_files/chr14.vg >chr14.gfa
/scratch/hllab/Juan/Latest_vg/vg convert -f vg_files/chr15.vg >chr15.gfa
/scratch/hllab/Juan/Latest_vg/vg convert -f vg_files/chr16.vg >chr16.gfa
/scratch/hllab/Juan/Latest_vg/vg convert -f vg_files/chr17.vg >chr17.gfa
/scratch/hllab/Juan/Latest_vg/vg convert -f vg_files/chr18.vg >chr18.gfa
/scratch/hllab/Juan/Latest_vg/vg convert -f vg_files/chr19.vg >chr19.gfa
/scratch/hllab/Juan/Latest_vg/vg convert -f vg_files/chr20.vg >chr20.gfa
/scratch/hllab/Juan/Latest_vg/vg convert -f vg_files/chr21.vg >chr21.gfa
/scratch/hllab/Juan/Latest_vg/vg convert -f vg_files/chr22.vg >chr22.gfa
/scratch/hllab/Juan/Latest_vg/vg convert -f vg_files/chrX.vg >chrX.gfa
/scratch/hllab/Juan/Latest_vg/vg convert -f vg_files/chrY.vg >chrY.gfa

exit
```
### Clean up
```bash
cd /scratch/hllab/Juan/Ixchel_Dev_Tests/Processs_hprc_v1_1_mc_chm13

# make .gfa directory
mkdir -p gfa
mv chr*.gfa gfa/
```

## Convert GFA to rGFA
### Link the gfa files to the working directory
```bash
cd /scratch/hllab/Juan/Ixchel_Dev_Tests/Processs_hprc_v1_1_mc_chm13

ln -s gfa/chr1.gfa .
ln -s gfa/chr10.gfa .
ln -s gfa/chr11.gfa .
ln -s gfa/chr12.gfa .
ln -s gfa/chr13.gfa .
ln -s gfa/chr14.gfa .
ln -s gfa/chr15.gfa .
ln -s gfa/chr16.gfa .
ln -s gfa/chr17.gfa .
ln -s gfa/chr18.gfa .
ln -s gfa/chr19.gfa .
ln -s gfa/chr20.gfa .
ln -s gfa/chr21.gfa .
ln -s gfa/chr22.gfa .
ln -s gfa/chr2.gfa .
ln -s gfa/chr3.gfa .
ln -s gfa/chr4.gfa .
ln -s gfa/chr5.gfa .
ln -s gfa/chr6.gfa .
ln -s gfa/chr7.gfa .
ln -s gfa/chr8.gfa .
ln -s gfa/chr9.gfa .
ln -s gfa/chrX.gfa .
ln -s gfa/chrY.gfa .
```
### Create params file
```bash
cd /scratch/hllab/Juan/Ixchel_Dev_Tests/Processs_hprc_v1_1_mc_chm13

paste \
<(find . -type l -name "chr*.gfa" | cut -d/ -f2) \
<(find . -type l -name "chr*.gfa" | cut -d/ -f2 | awk '{print "CHM13"}') \
<(find . -type l -name "chr*.gfa" | cut -d/ -f2 | sed s/'.gfa'/'.rgfa'/ ) \
> convert_gfa_to_rGFA_parameters.txt
```
### How many jobs?
```bash
cd /scratch/hllab/Juan/Ixchel_Dev_Tests/Processs_hprc_v1_1_mc_chm13

wc -l convert_gfa_to_rGFA_parameters.txt
```
```console
24 convert_gfa_to_rGFA_parameters.txt
```
### Launch jobs
```bash
cd /scratch/hllab/Juan/Ixchel_Dev_Tests/Processs_hprc_v1_1_mc_chm13

sbatch --array=1-24%24 --mem=12000 /scratch/hllab/Juan/JuanMacias_General_Code/Pangenomic/run_Batch_convert_GFA_to_rGFA.sh convert_gfa_to_rGFA_parameters.txt
```
```console
Submitted batch job 22931813
```
#### Check the status of the jobs
This is a step I do tailored to my compute environment; it is not essential for the process to work.
```bash
cd /scratch/hllab/Juan/Ixchel_Dev_Tests/Processs_hprc_v1_1_mc_chm13
bash /scratch/hllab/Juan/JuanMacias_General_Code/Job_Management/run_check_jobs_with_reportseff.sh 22931813 | cat
```
```console
Efficiency report for job ID: 22931813

        JobID    State       Elapsed  TimeEff   CPUEff   MemEff 
   22931813_1  COMPLETED    00:00:23   0.0%     65.2%     5.0%  
   22931813_2  COMPLETED    00:00:25   0.0%     64.0%     5.1%  
   22931813_3  COMPLETED    00:00:13   0.0%     61.5%     3.3%  
   22931813_4  COMPLETED    00:00:53   0.0%     67.9%    16.0%  
   22931813_5  COMPLETED    00:00:32   0.0%     65.6%     9.2%  
   22931813_6  COMPLETED    00:00:31   0.0%     64.5%     8.6%  
   22931813_7  COMPLETED    00:00:23   0.0%     65.2%     5.8%  
   22931813_8  COMPLETED    00:00:50   0.0%     68.0%    12.4%  
   22931813_9  COMPLETED    00:00:29   0.0%     62.1%     6.7%  
  22931813_10  COMPLETED    00:00:34   0.0%     67.6%     8.8%  
  22931813_11  COMPLETED    00:00:29   0.0%     65.5%     7.8%  
  22931813_12  COMPLETED    00:00:44   0.0%     65.9%    11.2%  
  22931813_13  COMPLETED    00:00:41   0.0%     65.9%     9.8%  
  22931813_14  COMPLETED    00:00:19   0.0%     63.2%     2.8%  
  22931813_15  COMPLETED    00:00:34   0.0%     67.6%    11.9%  
  22931813_16  COMPLETED    00:00:19   0.0%     63.2%     3.8%  
  22931813_17  COMPLETED    00:00:41   0.0%     65.9%    15.3%  
  22931813_18  COMPLETED    00:00:33   0.0%     66.7%    11.6%  
  22931813_19  COMPLETED    00:00:17   0.0%     64.7%     3.1%  
  22931813_20  COMPLETED    00:00:31   0.0%     67.7%    10.2%  
  22931813_21  COMPLETED    00:00:11   0.0%     54.5%     3.4%  
  22931813_22  COMPLETED    00:00:18   0.0%     66.7%     3.5%  
  22931813_23  COMPLETED    00:00:25   0.0%     64.0%     8.7%  
  22931813_24  COMPLETED    00:00:03   0.0%     33.3%     0.0%  

Log lengths:
   13 slurm-22931813_10.out
   13 slurm-22931813_11.out
   13 slurm-22931813_12.out
   13 slurm-22931813_13.out
   13 slurm-22931813_14.out
   13 slurm-22931813_15.out
   13 slurm-22931813_16.out
   13 slurm-22931813_17.out
   13 slurm-22931813_18.out
   13 slurm-22931813_19.out
   13 slurm-22931813_1.out
   13 slurm-22931813_20.out
   13 slurm-22931813_21.out
   13 slurm-22931813_22.out
   13 slurm-22931813_23.out
   13 slurm-22931813_24.out
   13 slurm-22931813_2.out
   13 slurm-22931813_3.out
   13 slurm-22931813_4.out
   13 slurm-22931813_5.out
   13 slurm-22931813_6.out
   13 slurm-22931813_7.out
   13 slurm-22931813_8.out
   13 slurm-22931813_9.out
  312 total
```
#### Cleanup
```bash
cd /scratch/hllab/Juan/Ixchel_Dev_Tests/Processs_hprc_v1_1_mc_chm13
mkdir -p Logs/convert_gfa_to_rGFA
mv slurm-*.out Logs/convert_gfa_to_rGFA/

# remove .gfa links
rm chr*.gfa

# move .rgfa files to rgfa directory
mkdir -p rgfa
mv chr*.rgfa rgfa/

# params file
mkdir -p Parameters
mv convert_gfa_to_rGFA_parameters.txt Parameters/

# create links to the .rgfa files
cd /scratch/hllab/Juan/Ixchel_Dev_Tests/Processs_hprc_v1_1_mc_chm13
ln -s rgfa/chr1.rgfa .
ln -s rgfa/chr10.rgfa .
ln -s rgfa/chr11.rgfa .
ln -s rgfa/chr12.rgfa .
ln -s rgfa/chr13.rgfa .
ln -s rgfa/chr14.rgfa .
ln -s rgfa/chr15.rgfa .
ln -s rgfa/chr16.rgfa .
ln -s rgfa/chr17.rgfa .
ln -s rgfa/chr18.rgfa .
ln -s rgfa/chr19.rgfa .
ln -s rgfa/chr20.rgfa .
ln -s rgfa/chr21.rgfa .
ln -s rgfa/chr22.rgfa .
ln -s rgfa/chr2.rgfa .
ln -s rgfa/chr3.rgfa .
ln -s rgfa/chr4.rgfa .
ln -s rgfa/chr5.rgfa .
ln -s rgfa/chr6.rgfa .
ln -s rgfa/chr7.rgfa .
ln -s rgfa/chr8.rgfa .
ln -s rgfa/chr9.rgfa .
ln -s rgfa/chrX.rgfa .
ln -s rgfa/chrY.rgfa .
```

## prepareGraphFiles
### Create params file
```bash
cd /scratch/hllab/Juan/Ixchel_Dev_Tests/Processs_hprc_v1_1_mc_chm13

paste \
<(find . -type l -name "chr*.rgfa" | cut -d/ -f2) \
<(find . -type l -name "chr*.rgfa" | cut -d/ -f2 | awk '{print "CHM13"}') \
> prepareGraphFiles_parameters.txt
```
### How many jobs?
```bash
cd /scratch/hllab/Juan/Ixchel_Dev_Tests/Processs_hprc_v1_1_mc_chm13

wc -l prepareGraphFiles_parameters.txt
```
```console
24 prepareGraphFiles_parameters.txt
```
### Launch jobs
```bash
cd /scratch/hllab/Juan/Ixchel_Dev_Tests/Processs_hprc_v1_1_mc_chm13
sbatch --array=1-24%24 --mem=12000 /scratch/hllab/Juan/JuanMacias_General_Code/Pangenomic/run_Batch_Ixchel_prepareGraphFiles.sh prepareGraphFiles_parameters.txt
```
```console
Submitted batch job 22931839
```
#### Check the status of the jobs
```bash
cd /scratch/hllab/Juan/Ixchel_Dev_Tests/Processs_hprc_v1_1_mc_chm13
bash /scratch/hllab/Juan/JuanMacias_General_Code/Job_Management/run_check_jobs_with_reportseff.sh 22931839 | cat
```
```console
Efficiency report for job ID: 22931839

        JobID    State       Elapsed  TimeEff   CPUEff   MemEff 
   22931839_1  COMPLETED    00:03:34   0.0%     91.1%    11.1%  
   22931839_2  COMPLETED    00:03:48   0.0%     90.4%    18.7%  
   22931839_3  COMPLETED    00:02:12   0.0%     92.4%     5.7%  
   22931839_4  COMPLETED    00:05:41   0.1%     93.0%    36.4%
   22931839_5  COMPLETED    00:05:39   0.1%     93.2%    34.7%
   22931839_6  COMPLETED    00:03:01   0.0%     92.3%    18.0%  
   22931839_7  COMPLETED    00:04:21   0.0%     92.0%    23.6%
   22931839_8  COMPLETED    00:04:38   0.0%     91.7%    20.7%
   22931839_9  COMPLETED    00:07:41   0.1%     92.4%    39.1%
  22931839_10  COMPLETED    00:00:53   0.0%     83.0%     1.3%  
  22931839_11  COMPLETED    00:02:36   0.0%     90.4%    11.1%  
  22931839_12  COMPLETED    00:07:44   0.1%     92.0%    47.0%
  22931839_13  COMPLETED    00:04:26   0.0%     91.7%    21.2%
  22931839_14  COMPLETED    00:03:08   0.0%     89.4%     6.3%  
  22931839_15  COMPLETED    00:02:59   0.0%     89.9%    21.7%
  22931839_16  COMPLETED    00:01:47   0.0%     89.7%     8.3%  
  22931839_17  COMPLETED    00:02:05   0.0%     89.6%    15.2%  
  22931839_18  COMPLETED    00:03:45   0.0%     91.1%    31.6%
  22931839_19  COMPLETED    00:03:44   0.0%     91.5%    30.4%
  22931839_20  COMPLETED    00:01:05   0.0%     89.2%     2.4%  
  22931839_21  COMPLETED    00:03:22   0.0%     91.1%    28.3%
  22931839_22  COMPLETED    00:03:59   0.0%     90.8%    30.5%
  22931839_23  COMPLETED    00:02:58   0.0%     91.0%    15.2%  
  22931839_24  COMPLETED    00:01:02   0.0%     87.1%     6.7%  

Log lengths:
     64 slurm-22931839_10.out
     64 slurm-22931839_11.out
     64 slurm-22931839_12.out
     64 slurm-22931839_13.out
     64 slurm-22931839_14.out
     64 slurm-22931839_15.out
     64 slurm-22931839_16.out
     64 slurm-22931839_17.out
     64 slurm-22931839_18.out
     64 slurm-22931839_19.out
     64 slurm-22931839_1.out
     64 slurm-22931839_20.out
     64 slurm-22931839_21.out
     64 slurm-22931839_22.out
     64 slurm-22931839_23.out
     64 slurm-22931839_24.out
     64 slurm-22931839_2.out
     64 slurm-22931839_3.out
     64 slurm-22931839_4.out
     64 slurm-22931839_5.out
     64 slurm-22931839_6.out
     64 slurm-22931839_7.out
     64 slurm-22931839_8.out
     64 slurm-22931839_9.out
   1536 total
```
#### Cleanup
```bash
cd /scratch/hllab/Juan/Ixchel_Dev_Tests/Processs_hprc_v1_1_mc_chm13
mkdir -p Logs/prepareGraphFiles
mv slurm-*.out Logs/prepareGraphFiles/
```

## convert split annotation files to precomputed conversion files
### Make parmeters file
Format: `{INPUTFILE} {SEGMENTSPICKLE} {QUERYSEGMENTSPICKLE} {LINKSPICKLE} {UpstreamLinksArrayFile} {DownstreamLinksArrayFile} {DoubleAnchorsFile}`
```bash
cd /scratch/hllab/Juan/Ixchel_Dev_Tests/Processs_hprc_v1_1_mc_chm13

paste \
<(find split_annotations/ -type f) \
<(find split_annotations/ -type f | while read -r line; do
    base=$(basename "$line")
    assembly=$(echo "$base" | cut -d '.' -f3 | cut -d '_' -f1)
    find . -type f -name "RefOnly.Segments.${assembly}.pkl" -exec basename {} \;
done) \
<(find split_annotations/ -type f | while read -r line; do
    base=$(basename "$line")
    assembly=$(echo "$base" | cut -d '.' -f3 | cut -d '_' -f1)
    find . -type f -name "QueryOnly.Segments.${assembly}.pkl" -exec basename {} \;
done) \
<(find split_annotations/ -type f | while read -r line; do
    base=$(basename "$line")
    assembly=$(echo "$base" | cut -d '.' -f3 | cut -d '_' -f1)
    find . -type f -name "FilteredLinks.Links.${assembly}.pkl" -exec basename {} \;
done) \
<(find split_annotations/ -type f | while read -r line; do
    base=$(basename "$line")
    assembly=$(echo "$base" | cut -d '.' -f3 | cut -d '_' -f1)
    find . -type f -name "UpstreamArray.RefOnly.Segments.${assembly}.pkl" -exec basename {} \;
done) \
<(find split_annotations/ -type f | while read -r line; do
    base=$(basename "$line")
    assembly=$(echo "$base" | cut -d '.' -f3 | cut -d '_' -f1)
    find . -type f -name "DownstreamArray.RefOnly.Segments.${assembly}.pkl" -exec basename {} \;
done) \
<(find split_annotations/ -type f | while read -r line; do
    base=$(basename "$line")
    assembly=$(echo "$base" | cut -d '.' -f3 | cut -d '_' -f1)
    find . -type f -name "DoubleAnchored.FilteredLinks.Links.${assembly}.pkl" -exec basename {} \;
done) \
> Convert_split_annotations_parameters.txt
```
### How many jobs?
```bash
wc -l Convert_split_annotations_parameters.txt
```
```console
1271 Convert_split_annotations_parameters.txt
```
### Run the jobs
```bash
cd /scratch/hllab/Juan/Ixchel_Dev_Tests/Processs_hprc_v1_1_mc_chm13
sbatch  --array=1-1271%175 --mem=12G /scratch/hllab/Juan/Ixchel/SourceCode/bash_Scripts_For_Stepwise_Processing/run_Batch_PreConvertAnnotations.sh Convert_split_annotations_parameters.txt
```
```console
Submitted batch job 22931867
```
### Check the status of the jobs
```bash
cd /scratch/hllab/Juan/Ixchel_Dev_Tests/Processs_hprc_v1_1_mc_chm13
bash /scratch/hllab/Juan/JuanMacias_General_Code/Job_Management/run_check_jobs_with_reportseff.sh -s 22931867
```
```console
Efficiency report for job ID: 22931867

States:
  COMPLETED: 1271

Metric          min     max     mean    stddev
Elapsed(min)    0.17    486.08  6.96    23.49
TimeEff(%)      0.0     4.8     0.06    0.24
CPUEff(%)       78.9    99.8    98.40   1.71
MemEff(%)       0.0     8.7     32.60   12.02

Log-file line-count frequencies (lines → #files):
    19 → 1271
```
### Cleanup
```bash
cd /scratch/hllab/Juan/Ixchel_Dev_Tests/Processs_hprc_v1_1_mc_chm13
mkdir -p Logs/Convert_split_annotations
mv slurm-* Logs/Convert_split_annotations/

mkdir -p split_conversions
mv split_annotations/*.converted split_conversions/

mkdir -p Parameters
mv Convert_split_annotations_parameters.txt Parameters/
mv prepareGraphFiles_parameters.txt Parameters/

# remove sym links to the .rgfa files
rm chr*.rgfa

# move preprocessing intermediate files to a separate directory
mkdir -p Preprocessing_Intermediate_Files
mv Annotations.* Preprocessing_Intermediate_Files/
mv DoubleAnchored.*pkl Preprocessing_Intermediate_Files/
mv DownstreamArray.*pkl Preprocessing_Intermediate_Files/
mv FilteredLinks.*pkl Preprocessing_Intermediate_Files/
mv FilteredLinks.*rgfa Preprocessing_Intermediate_Files/
mv Links.*.rgfa Preprocessing_Intermediate_Files/
mv QueryOnly.*pkl Preprocessing_Intermediate_Files/
mv QueryOnly.*rgfa Preprocessing_Intermediate_Files/
mv RefOnly.*pkl Preprocessing_Intermediate_Files/
mv RefOnly.*rgfa Preprocessing_Intermediate_Files/
mv Segments.*rgfa Preprocessing_Intermediate_Files/
mv UpstreamArray.*pkl Preprocessing_Intermediate_Files/
```
#### Tar and compress the `Preprocessing_Intermediate_Files`
```bash
cd /scratch/hllab/Juan/Ixchel_Dev_Tests/Processs_hprc_v1_1_mc_chm13

srun --mem=16000 --cpus-per-task=20 -J interactive -p interactive --pty /bin/bash -l
tar -cvzf Preprocessing_Intermediate_Files.tar.gz Preprocessing_Intermediate_Files/
exit
```
```console
...
Preprocessing_Intermediate_Files/Links.chr5.rgfa
Preprocessing_Intermediate_Files/Segments.chr11.rgfa
Preprocessing_Intermediate_Files/Segments.chr8.rgfa
Preprocessing_Intermediate_Files/QueryOnly.Segments.chr10.rgfa
Preprocessing_Intermediate_Files/FilteredLinks.Links.chr1.rgfa
Preprocessing_Intermediate_Files/FilteredLinks.Links.chr8.pkl
Preprocessing_Intermediate_Files/RefOnly.Segments.chr10.pkl
Preprocessing_Intermediate_Files/FilteredLinks.Links.chr12.pkl
Preprocessing_Intermediate_Files/UpstreamArray.RefOnly.Segments.chr10.pkl
Preprocessing_Intermediate_Files/Links.chr2.rgfa
Preprocessing_Intermediate_Files/QueryOnly.Segments.chr11.pkl
```

## prepare lookups
### Cat the split_conversions
```bash
cd /scratch/hllab/Juan/Ixchel_Dev_Tests/Processs_hprc_v1_1_mc_chm13

srun --mem=16000 --cpus-per-task=1 -J interactive -p interactive --pty /bin/bash -l
cat split_conversions/*.converted >Annotations.converted
exit
```
### Breakdown of flags
```bash
cd /scratch/hllab/Juan/Ixchel_Dev_Tests/Processs_hprc_v1_1_mc_chm13

srun --mem=32000 --cpus-per-task=1 -J interactive -p interactive --pty /bin/bash -l

cut -f8 Annotations.converted | sort | uniq -c

exit
```
```console
 466678 1
2419485 100
  39737 102
 639097 114
 193120 118
1156870976 3
9421558 34
 548840 38
 242873 42
9632499 46
  17107 5
13586989 7
58918028 96
5036482 98
```
| %      | Count          | Code |
|--------|----------------|------|
| 0.0%   | 466,678        | 1    |
| 0.2%   | 2,419,485      | 100  |
| 0.0%   | 39,737         | 102  |
| 0.1%   | 639,097        | 114  |
| 0.0%   | 193,120        | 118  |
| 92.0%  | 1,156,870,976  | 3    |
| 0.7%   | 9,421,558      | 34   |
| 0.0%   | 548,840        | 38   |
| 0.0%   | 242,873        | 42   |
| 0.8%   | 9,632,499      | 46   |
| 0.0%   | 17,107         | 5    |
| 1.1%   | 13,586,989     | 7    |
| 4.7%   | 58,918,028     | 96   |
| 0.4%   | 5,036,482      | 98   |
### Build the database
```bash
cd /scratch/hllab/Juan/Ixchel_Dev_Tests/Processs_hprc_v1_1_mc_chm13

sbatch /scratch/hllab/Juan/JuanMacias_General_Code/Pangenomic/run_Ixchel_Build_DB.sh Annotations.converted Annotations.converted.db
```
```console
Submitted batch job 22933394
```
#### Check the status of the jobs
```bash
cd /scratch/hllab/Juan/Ixchel_Dev_Tests/Processs_hprc_v1_1_mc_chm13
bash /scratch/hllab/Juan/JuanMacias_General_Code/Job_Management/run_check_jobs_with_reportseff.sh 22933394
```
```console
Efficiency report for job ID: 22933394

     JobID    State       Elapsed  TimeEff   CPUEff   MemEff 
  22933394  COMPLETED    02:38:29   1.6%     65.7%    31.8%

Log lengths:
258 slurm-22933394.out
```
### Cleanup
```bash
cd /scratch/hllab/Juan/Ixchel_Dev_Tests/Processs_hprc_v1_1_mc_chm13
rm -r Preprocessing_Intermediate_Files

mkdir -p Logs/Build_DB
mv slurm-*.out Logs/Build_DB/

# compress `Annotations.converted`
srun --mem=16000 --cpus-per-task=8 -J interactive -p interactive --pty /bin/bash -l
pigz Annotations.converted
exit

cp Annotations.converted.db HPRC_v1_1_mc_CHM14_Annotations.converted.db

# compress database that will be made available
srun --mem=16000 --cpus-per-task=8 -J interactive -p interactive --pty /bin/bash -l
pigz HPRC_v1_1_mc_CHM14_Annotations.converted.db
exit
```

The `Annotations.converted.db` is the main file needed for the conversion process.
