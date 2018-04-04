#!/bin/bash

name="example"               #-- any name
outDir="/example/outFolder"  #-- any dir
bed="/example/target.bed"           #-- any bed file

hg19="/gscmnt/gc2737/ding/hsun/data/GRCh37-lite/GRCh37-lite.fa"
mm10="/gscmnt/gc2737/ding/hsun/data/ensemble_v91/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa"
chain10="/gscmnt/gc2737/ding/hsun/data/liftOver/hg19ToMm10.over.chain"

mkdir -p $outDir/log

gmt somatic filter-mouse-bases --chain-file=$chain10 --human-reference=$hg19 --mouse-reference=$mm10 --variant-file=$bed --filtered-file=$outDir/log/$name.mouse.hg19toMm10.out --output-file=$outDir/$name.hg19toMm10.permissive.out --permissive > $outDir/log/$name.mouseFilter.hg19toMm10.log


## target.bed
#1	121009	121009	C	T
#1	234308	234308	A	G

