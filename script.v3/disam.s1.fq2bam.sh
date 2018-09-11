#!/bin/bash

# 2018/09/08
# Hua Sun

# sh disam.s1.fq2bam.sh -C contig.ini -N test -S human -1 name.fq1.gz -2 name.fq2.gz -O /path/outDir
# memory 10 Gb

# getOptions
while getopts "C:N:S:1:2:O:" opt; do
  case $opt in
    C)
      CONFIG=$OPTARG
      ;;
    N)
      NAME=$OPTARG
      ;;
    S)
      SPECIES=$OPTARG
      ;;
    1)
      FQ1=$OPTARG
      ;;
    2)
      FQ2=$OPTARG
      ;;    
    O)
      OUTDIR=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done


source $CONFIG

if [ ! -d $OUTDIR ]; then
	echo "[ERROR] The $OUTDIR not exists!" >&2
	exit 1
fi

if [ -z "$NAME" ]; then
	echo "[ERROR] The Name is empty!" >&2
	exit 1
fi

GENOME=''

# set species
if [ "$SPECIES" = "human" ]; then
	GENOME=$REF_HUMAN
elif [ "$SPECIES" = "mouse" ]; then
	GENOME=$REF_MOUSE
else
	echo "[ERROR] The -S should be human or mouse!" >&2
	exit 1
fi

if [ ! -f $GENOME ]; then
	echo "[ERROR] The $GENOME not exists!" >&2
	exit 1
fi

# check file
if [ -f $OUTDIR/$NAME.$SPECIES.bam ]; then
    echo "The $OUTDIR/$NAME.$SPECIES.bam file exists!" >&2
    exit 1
fi

# bwa hg/mm
$BWA mem -t 8 -M -R "@RG\tID:$NAME\tSM:$NAME\tPL:illumina\tLB:$NAME.lib\tPU:$NAME.unit" $GENOME $FQ1 $FQ2 > $OUTDIR/$NAME.$SPECIES.sam
$SAMTOOLS view -Sbh $OUTDIR/$NAME.$SPECIES.sam > $OUTDIR/$NAME.$SPECIES.bam

# sort bam by natural name
$SAMTOOLS sort -m 1G -@ 6 -o $OUTDIR/$NAME.$SPECIES.sort.bam -n $OUTDIR/$NAME.$SPECIES.bam

