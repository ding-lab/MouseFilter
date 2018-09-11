#!/bin/bash

# 2018/09/08
# Hua Sun

# sh disam.s2.newBam.sh -C contig.ini -N test -O /path/inputDir
# memory 18 Gb
# the -O should be same with 'disam.s1.fq2bam.sh' folder

# getOptions
while getopts "C:N:O:" opt; do
  case $opt in
    C)
      CONFIG=$OPTARG
      ;;
    N)
      NAME=$OPTARG
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

if [ ! -f $OUTDIR/$NAME.human.sort.bam ]; then
	echo "[ERROR] No $OUTDIR/$NAME.human.sort.bam file!" >&2
	exit 1
fi

if [ -z "$NAME" ]; then
	echo "[ERROR] The Name is empty!" >&2
	exit 1
fi


# Disambiguate (mouse-filter)
$DISAMBIGUATE -s $NAME -o $OUTDIR -a bwa $OUTDIR/$NAME.human.sort.bam $OUTDIR/$NAME.mouse.sort.bam

# re-create fq
$SAMTOOLS sort -m 1G -@ 6 -o $OUTDIR/$NAME.disam.sortbyname.bam -n $OUTDIR/$NAME.disambiguatedSpeciesA.bam
$SAMTOOLS fastq $OUTDIR/$NAME.disam.sortbyname.bam -1 $OUTDIR/$NAME.disam_1.fastq.gz -2 $OUTDIR/$NAME.disam_2.fastq.gz


# mapping to human reference and create to bam
$BWA mem -t 8 -M -R "@RG\tID:$NAME\tSM:$NAME\tPL:illumina\tLB:$NAME.lib\tPU:$NAME.unit" $REF_HUMAN $OUTDIR/$NAME.disam_1.fastq.gz $OUTDIR/$NAME.disam_2.fastq.gz > $OUTDIR/$NAME.disam.reAlign.sam

# sort
$JAVA -Xmx16G -jar $PICARD SortSam \
   I=$OUTDIR/$NAME.disam.reAlign.sam \
   O=$OUTDIR/$NAME.disam.reAlign.bam \
   SORT_ORDER=coordinate

# remove-duplication
$JAVA -Xmx16G -jar $PICARD MarkDuplicates \
   I=$OUTDIR/$NAME.disam.reAlign.bam \
   O=$OUTDIR/$NAME.disam.reAlign.remDup.bam \
   REMOVE_DUPLICATES=true \
   M=$OUTDIR/$NAME.disam.reAlign.remDup.metrics.txt

# index bam
$SAMTOOLS index $OUTDIR/$NAME.disam.reAlign.remDup.bam

