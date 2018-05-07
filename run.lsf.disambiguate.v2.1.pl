=head1 Date
  4/30/2018 updated
  Hua Sun

=head1 DESCRIPTION

	A testing version

=head1 USAGE

  perl run.lsf.disambiguate.pl -t dna -q long -f bam.table -outDir /path/data
  perl run.lsf.disambiguate.pl -t dna -q long -f bam.table -outDir /path/data -run
  perl run.lsf.disambiguate.pl -t rna -q long -f bam.table -outDir /path/data -run

  NOTE: #bam.table
        name	sampleID /path/*.bam

=head1 OPTIONS
  
  -t        [str]  #dna/rna
  -q        [str]  #bqueues (long)
  -f        [file] #bamList file
  -outDir   [str]  #result path (in the path including several detail sample folders)
  -run             #submit job (bsub)
    
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

my ($name,$outFolder);
foreach my $str (@data){
  chomp($str);
  next if ($str eq '');
  
  ($name,$outFolder) = &Create_DNA_Script($que, $str, $outDir) if ($type eq 'dna');
  ($name,$outFolder) = &Create_RNA_Script($que, $str, $outDir) if ($type eq 'rna');
  
  print "Create Job: $name\n";
  
  system("bsub < $outFolder/lsf.$name.script") if (defined $run);
  
  sleep 1;
}

exit;

##########################################################

##### DNA mouse-filter
sub Create_DNA_Script
{
  my ($que, $str, $outDir)=@_;
  my ($name, $sample, $bam) = split("\t", $str); 
  
  my $outFolder = "$outDir/$name";
  system("mkdir -p $outFolder");  #create new directory
  
  open my ($OUT), ">$outFolder/lsf.$name.script";
  
  print $OUT <<EOF
#!/bin/bash
#BSUB -q $que
#BSUB -R 'select[mem>16000] rusage[mem=16000] span[hosts=1]'
#BSUB -M 16000000
#BSUB -n 1
#BSUB -J "$name"
#BSUB -g "/hsun"
#BSUB -eo $outFolder/$name.err
#BSUB -oo $outFolder/$name.log

JAVA="/gscmnt/gc2737/ding/hsun/software/jre1.8.0_152/bin/java"
picard="/gscmnt/gc2737/ding/hsun/software/picard.jar"
bwa="/gscmnt/gc2737/ding/hsun/software/bwa-0.7.17/bwa"
samtools="/gscmnt/gc2737/ding/hsun/software/miniconda2/bin/samtools"
disambiguate="/gscmnt/gc2737/ding/hsun/software/miniconda2/bin/ngs_disambiguate"

hg_genome="/gscmnt/gc2737/ding/hsun/data/GRCh37-lite/GRCh37-lite.fa"
mm_genome="/gscmnt/gc2737/ding/hsun/data/ensemble_v91/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa"

# bam2fq
\$samtools sort -m 1G -@ 6 -o $outFolder/$name.sortbyname.bam -n $bam
\$samtools fastq $outFolder/$name.sortbyname.bam -1 $outFolder/$name\_1.fastq.gz -2 $outFolder/$name\_2.fastq.gz

# bwa hg/mm
\$bwa mem -t 8 -M -R "\@RG\\tID:$sample\\tSM:$sample\\tPL:illumina\\tLB:$sample.lib\\tPU:$sample.unit" \$hg_genome $outFolder/$name\_1.fastq.gz $outFolder/$name\_2.fastq.gz > $outFolder/human.sam
\$samtools view -Sbh $outFolder/human.sam > $outFolder/human.bam
\$bwa mem -t 8 -M -R "\@RG\\tID:$sample\\tSM:$sample\\tPL:illumina\\tLB:$sample.lib\\tPU:$sample.unit" \$mm_genome $outFolder/$name\_1.fastq.gz $outFolder/$name\_2.fastq.gz > $outFolder/mouse.sam
\$samtools view -Sbh $outFolder/mouse.sam > $outFolder/mouse.bam

# sort bam by natural name
\$samtools sort -m 1G -@ 6 -o $outFolder/human.sort.bam -n $outFolder/human.bam
\$samtools sort -m 1G -@ 6 -o $outFolder/mouse.sort.bam -n $outFolder/mouse.bam

# Disambiguate (mouse-filter)
\$disambiguate -s $name -o $outFolder -a bwa $outFolder/human.sort.bam $outFolder/mouse.sort.bam

# re-create fq
\$samtools sort -m 1G -@ 6 -o $outFolder/$name.disam.sortbyname.bam -n $outFolder/$name.disambiguatedSpeciesA.bam
\$samtools fastq $outFolder/$name.disam.sortbyname.bam -1 $outFolder/$name.disam_1.fastq.gz -2 $outFolder/$name.disam_2.fastq.gz

# bwa hg
\$bwa mem -t 8 -M -R "\@RG\\tID:$sample\\tSM:$sample\\tPL:illumina\\tLB:$sample.lib\\tPU:$sample.unit" \$hg_genome $outFolder/$name.disam_1.fastq.gz $outFolder/$name.disam_2.fastq.gz > $outFolder/$name.disam.reAligned.sam

# sort
\$JAVA -Xmx8G -jar \$picard SortSam \\
   I=$outFolder/$name.disam.reAligned.sam \\
   O=$outFolder/$name.disam.reAligned.bam \\
   SORT_ORDER=coordinate

# remove-duplication
\$JAVA -Xmx8G -jar \$picard MarkDuplicates \\
   I=$outFolder/$name.disam.reAligned.bam \\
   O=$outFolder/$name.disam.reAligned.remDup.bam \\
   REMOVE_DUPLICATES=true \\
   M=$outFolder/picard.disam.reAligned.remdup.metrics.txt

# index bam
\$samtools index $outFolder/$name.disam.reAligned.remDup.bam

EOF
    ;

  return($name,$outFolder);
}



##### RNA mouse-filter
sub Create_RNA_Script
{
  my ($que, $str, $outDir)=@_;
  my ($name, $sample, $bam) = split("\t", $str);  # $sample does not use in RNA filter code
  
  my $outFolder = "$outDir/$name";
  system("mkdir -p $outFolder");  #create new directory
  
  open my ($OUT), ">$outFolder/lsf.$name.script";
  
  print $OUT <<EOF
#!/bin/bash
#BSUB -q $que
#BSUB -R 'select[mem>64000] rusage[mem=64000] span[hosts=1]'
#BSUB -M 64000000
#BSUB -n 1
#BSUB -J "$name"
#BSUB -g "/hsun"
#BSUB -eo $outFolder/$name.err
#BSUB -oo $outFolder/$name.log

JAVA="/gscmnt/gc2737/ding/hsun/software/jre1.8.0_152/bin/java"
picard="/gscmnt/gc2737/ding/hsun/software/picard.jar"
star="/gscmnt/gc2737/ding/hsun/software/miniconda2/bin/STAR"
samtools="/gscmnt/gc2737/ding/hsun/software/miniconda2/bin/samtools"
disambiguate="/gscmnt/gc2737/ding/hsun/software/miniconda2/bin/ngs_disambiguate"

hg_genomeDir="/gscmnt/gc2737/ding/hsun/data/ensemble_v91/GRCh37_star_genomeDir"
hg_gtf="/gscmnt/gc2737/ding/hsun/data/ensemble_v91/Homo_sapiens.GRCh37.87.gtf"
hg_genome="/gscmnt/gc2737/ding/hsun/data/ensemble_v91/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa"

mm_genomeDir="/gscmnt/gc2737/ding/hsun/data/ensemble_v91/mm10_star_genomeDir"
mm_gtf="/gscmnt/gc2737/ding/hsun/data/ensemble_v91/Mus_musculus.GRCm38.91.gtf"
mm_genome="/gscmnt/gc2737/ding/hsun/data/ensemble_v91/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa"


## bam2fq
\$samtools sort -m 1G -@ 6 -o $outFolder/$name.sortbyname.bam -n $bam
\$samtools fastq $outFolder/$name.sortbyname.bam -1 $outFolder/$name\_1.fastq.gz -2 $outFolder/$name\_2.fastq.gz

# align sequencing reads to the genome
\$star --runThreadN 12 --genomeDir \$hg_genomeDir --sjdbGTFfile \$hg_gtf --sjdbOverhang 100 --readFilesIn $outFolder/$name\_1.fastq.gz $outFolder/$name\_2.fastq.gz --outFileNamePrefix $outFolder/human. --outSAMtype BAM Unsorted --twopassMode Basic --outSAMattributes All --genomeLoad NoSharedMemory --readFilesCommand zcat
\$star --runThreadN 12 --genomeDir \$mm_genomeDir --sjdbGTFfile \$mm_gtf --sjdbOverhang 100 --readFilesIn $outFolder/$name\_1.fastq.gz $outFolder/$name\_2.fastq.gz --outFileNamePrefix $outFolder/mouse. --outSAMtype BAM Unsorted --twopassMode Basic --outSAMattributes All --genomeLoad NoSharedMemory --readFilesCommand zcat

# sort bam by natural name
\$samtools sort -m 1G -@ 6 -o $outFolder/human.sort.bam -n $outFolder/human.Aligned.out.bam
\$samtools sort -m 1G -@ 6 -o $outFolder/mouse.sort.bam -n $outFolder/mouse.Aligned.out.bam

# Disambiguate (mouse-filter)
\$disambiguate -s $name -o $outFolder -a star $outFolder/human.sort.bam $outFolder/mouse.sort.bam

# re-create fq
\$samtools sort -m 1G -@ 6 -o $outFolder/$name.disam.sortbyname.bam -n $outFolder/$name.disambiguatedSpeciesA.bam
\$samtools fastq $outFolder/$name.disam.sortbyname.bam -1 $outFolder/$name.disambiguated.1.fastq.gz -2 $outFolder/$name.disambiguated.2.fastq.gz


EOF
    ;

  return($name,$outFolder);
}

