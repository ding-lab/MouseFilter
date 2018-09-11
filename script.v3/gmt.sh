#!/bin/bash
# 2018/09/09
# Hua Sun

# sh gmt.sh -C contig.ini -N sampleName -B sampleName.bed -O /path/outDir
# memory 1 Gb

# getOptions
while getopts "C:N:B:O:" opt; do
  case $opt in
    C)
      CONFIG=$OPTARG
      ;;
    N)
      NAME=$OPTARG
      ;;
    B)
      BED=$OPTARG
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

if [ ! -f $BED ]; then
	echo "[ERROR] The $BED not exists!" >&2
	exit 1
fi

mkdir -p $OUTDIR/log

gmt somatic filter-mouse-bases --chain-file=$LIFTOVER_CHAIN --human-reference=$REF_HUMAN --mouse-reference=$REF_MOUSE --variant-file=$BED --filtered-file=$OUTDIR/log/$NAME.mouse.hg19toMm10.out --output-file=$OUTDIR/$NAME.hg19toMm10.permissive.out --permissive > $OUTDIR/log/$NAME.mouseFilter.hg19toMm10.log


## bed file
#1	121009	121009	C	T
#1	234308	234308	A	G

