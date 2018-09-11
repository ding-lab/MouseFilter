#!/bin/bash

# 2018/09/08
# Hua Sun

# sh disambiguate.full-step.from_fq.sh -C contig.ini -N test -1 name.fq1.gz -2 name.fq2.gz -O /path/outDir
# memory 18 Gb

# getOptions
while getopts "C:N:1:2:O:" opt; do
  case $opt in
    C)
      CONFIG=$OPTARG
      ;;
    N)
      NAME=$OPTARG
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


# human
# bwa 
$BWA mem -t 8 -M -R "@RG\tID:$NAME\tSM:$NAME\tPL:illumina\tLB:$NAME.lib\tPU:$NAME.unit" $REF_HUMAN $FQ1 $FQ2 > $OUTDIR/$NAME.human.sam
$SAMTOOLS view -Sbh $OUTDIR/$NAME.human.sam > $OUTDIR/$NAME.human.bam

# sort bam by natural name
$SAMTOOLS sort -m 1G -@ 6 -o $OUTDIR/$NAME.human.sort.bam -n $OUTDIR/$NAME.human.bam


# mouse
$BWA mem -t 8 -M -R "@RG\tID:$NAME\tSM:$NAME\tPL:illumina\tLB:$NAME.lib\tPU:$NAME.unit" $REF_MOUSE $FQ1 $FQ2 > $OUTDIR/$NAME.mouse.sam
$SAMTOOLS view -Sbh $OUTDIR/$NAME.mouse.sam > $OUTDIR/$NAME.mouse.bam

# sort bam by natural name
$SAMTOOLS sort -m 1G -@ 6 -o $OUTDIR/$NAME.mouse.sort.bam -n $OUTDIR/$NAME.mouse.bam


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

