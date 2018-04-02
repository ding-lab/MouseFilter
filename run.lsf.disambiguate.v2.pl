=head1 Date
  3/9/2018 updated
  Hua Sun

=head1 USAGE

  perl run.lsf.disambiguate.pl -t dna -q long -f bam.table -outDir /path/data
  perl run.lsf.disambiguate.pl -t dna -q long -f bam.table -outDir /path/data -run

  NOTE: #bam.table
        key	sampleName /path/*.bam

=head1 OPTIONS
  
  -t        [str]  #dna/rna
  -q        [str]  #bqueues (long / ding-lab)
  -f        [file] #bamList file
  -outDir   [str]  #result path (in the path including several detail sample folders)
  -run             #submit job (bsub)

=head1 INPUT
  key	/path/*.bam
    
=cut

use strict;
use Getopt::Long;

my $type = 'dna';
my $que = 'long';
my $fileList = '';
my $outDir = '.';
my $run;
GetOptions(
  "t:s" => \$type,
  "q:s" => \$que,
  "f:s" => \$fileList,
  "outDir:s" => \$outDir,
  "run" => \$run
);

die `pod2text $0` if ($outDir eq '' && $fileList eq '');

my @data = readpipe("cat $fileList");

system("mkdir -p $outDir");

my ($key,$outFolder);
foreach my $str (@data){
  chomp($str);
  next if ($str eq '');
  
  ($key,$outFolder) = &Create_DNA_Script($que, $str, $outDir) if ($type eq 'dna');
  
  print "Create Job: $key\n";
  
  system("bsub < $outFolder/lsf.$key.script") if (defined $run);
  
  sleep 1;
}

exit;

##########################################################

sub Create_DNA_Script
{
  my ($que, $str, $outDir)=@_;
  my ($key, $sample, $bam) = split("\t", $str); 
  
  my $outFolder = "$outDir/$key";
  system("mkdir -p $outFolder");  #create new directory
  
  open my ($OUT), ">$outFolder/lsf.$key.script";
  
  print $OUT <<EOF
#!/bin/bash
#BSUB -q $que
#BSUB -R 'select[mem>16000] rusage[mem=16000] span[hosts=1]'
#BSUB -M 16000000
#BSUB -n 1
#BSUB -J "$key"
#BSUB -g "/hsun"
#BSUB -eo $outFolder/$key.err
#BSUB -oo $outFolder/$key.log

JAVA="/gscmnt/gc2737/ding/hsun/software/jre1.8.0_152/bin/java"
picard="/gscmnt/gc2737/ding/hsun/software/picard.jar"
bwa="/gscmnt/gc2737/ding/hsun/software/bwa-0.7.17/bwa"
samtools="/gscmnt/gc2737/ding/hsun/software/miniconda2/bin/samtools"
disambiguate="/gscmnt/gc2737/ding/hsun/software/miniconda2/bin/ngs_disambiguate"

hg_genome="/gscmnt/gc2737/ding/hsun/data/GRCh37-lite/GRCh37-lite.fa"
mm_genome="/gscmnt/gc2737/ding/hsun/data/ensemble_v91/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa"

# bam2fq
\$samtools sort -m 1G -@ 6 -o $outFolder/$key.sortbyname.bam -n $bam
\$samtools fastq $outFolder/$key.sortbyname.bam -1 $outFolder/$key\_1.fastq.gz -2 $outFolder/$key\_2.fastq.gz

# bwa hg/mm
\$bwa mem -t 8 -M -R "\@RG\\tID:$sample\\tSM:$sample\\tPL:illumina\\tLB:$sample.lib\\tPU:$sample.unit" \$hg_genome $outFolder/$key\_1.fastq.gz $outFolder/$key\_2.fastq.gz > $outFolder/human.sam
\$samtools view -Sbh $outFolder/human.sam > $outFolder/human.bam
\$bwa mem -t 8 -M -R "\@RG\\tID:$sample\\tSM:$sample\\tPL:illumina\\tLB:$sample.lib\\tPU:$sample.unit" \$mm_genome $outFolder/$key\_1.fastq.gz $outFolder/$key\_2.fastq.gz > $outFolder/mouse.sam
\$samtools view -Sbh $outFolder/mouse.sam > $outFolder/mouse.bam

# sort bam by natural name
\$samtools sort -m 1G -@ 6 -o $outFolder/human.sort.bam -n $outFolder/human.bam
\$samtools sort -m 1G -@ 6 -o $outFolder/mouse.sort.bam -n $outFolder/mouse.bam

# Disambiguate (mouse-filter)
\$disambiguate -s $key -o $outFolder -a bwa $outFolder/human.sort.bam $outFolder/mouse.sort.bam

# re-create fq
\$samtools sort -m 1G -@ 6 -o $outFolder/$key.disam.sortbyname.bam -n $outFolder/$key.disambiguatedSpeciesA.bam
\$samtools fastq $outFolder/$key.disam.sortbyname.bam -1 $outFolder/$key.disam_1.fastq.gz -2 $outFolder/$key.disam_2.fastq.gz

# bwa hg
\$bwa mem -t 8 -M -R "\@RG\\tID:$sample\\tSM:$sample\\tPL:illumina\\tLB:$sample.lib\\tPU:$sample.unit" \$hg_genome $outFolder/$key.disam_1.fastq.gz $outFolder/$key.disam_2.fastq.gz > $outFolder/$key.disam.reAligned.sam

# sort
\$JAVA -Xmx8G -jar \$picard SortSam \\
   I=$outFolder/$key.disam.reAligned.sam \\
   O=$outFolder/$key.disam.reAligned.bam \\
   SORT_ORDER=coordinate

# remove-duplication
\$JAVA -Xmx8G -jar \$picard MarkDuplicates \\
   I=$outFolder/$key.disam.reAligned.bam \\
   O=$outFolder/$key.disam.reAligned.remDup.bam \\
   REMOVE_DUPLICATES=true \\
   M=$outFolder/picard.disam.reAligned.remdup.metrics.txt

# index bam
\$samtools index $outFolder/$key.disam.reAligned.remDup.bam

EOF
    ;

  return($key,$outFolder);
}


