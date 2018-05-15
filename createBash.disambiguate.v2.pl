=head1 Date
  
  4/30/2018; 3/29/2018
  Hua Sun

=head1 DESCRIPTION
	The script only run MGI-LSF system (-run).
	If let it run in normal Linux, it can manually run the script (Do not use -run).
	
	Disambiguate mouse reads

=head1 USAGE

	perl createBash.disambiguate.v2.pl -f folderList -p bam2fq
  perl createBash.disambiguate.v2.pl -f folderList -p fq-sortBam -s human
  perl createBash.disambiguate.v2.pl -f folderList -p fq-sortBam -s mouse
	
	perl createBash.disambiguate.v2.pl -f folderList -p bam2fq -run
	perl createBash.disambiguate.v2.pl -f folderList -p fq2disam -a bwa -t dna -run
	perl createBash.disambiguate.v2.pl -f folderList -p fq2disam -a star -t rna -run
	
  NOTE: //folderList (240_wxs folder including .bam or .fastq.gz)
	e.g:
	/data/240_wxs
	/data/250_wxs

=head1 OPTIONS
  
  -t  [str]  # dna(default)/rna
  -a  [str]  # algorithm bwa(default)
  -s  [str]  # species human(default) | mouse
  -f  [file] # folderList
  -run       # submit job
  -p  [str]  #program-code 
			        bam2fq
			        fq-sortBam
			        disam
			        newBam        #create New bam after Disambiguate using *.disambiguatedSpeciesA.bam
			        fq2disam ( -t dna -a bwa )
			        fq2disam ( -t rna -a star )
			        
			        
=head1 folderList
	e.g: (240_wxs folder including same folder name of 240_wxs.bam or 240_wxs.fastq.gz)
	/data/240_wxs
	/data/250_wxs
	
    
=cut

use strict;
use Cwd qw(cwd);
use Getopt::Long;

my $type = 'dna';
my $algorithm = 'bwa';
my $species = 'human';
my $folderList = '';
my $outDir = '.';
my $program = '';
my $run;
GetOptions(
  "t:s" => \$type,
  "s:s" => \$species,
  "a:s" => \$algorithm,
  "f:s" => \$folderList,
  "outDir:s" => \$outDir,
  "p:s" => \$program,
  "run" => \$run
);

die `pod2text $0!` if ($folderList eq '');

##-- COLOR CODE
my $red = "\e[31m";
my $gray = "\e[37m";
my $yellow = "\e[33m";
my $green = "\e[32m";
my $purple = "\e[35m";
my $cyan = "\e[36m";
my $normal = "\e[0m";


##-- Setting
our $JAVA = "/gscmnt/gc2737/ding/hsun/software/jre1.8.0_152/bin/java";
our $PICARD = "/gscmnt/gc2737/ding/hsun/software/picard.jar";
our $BWA = "/gscmnt/gc2737/ding/hsun/software/bwa-0.7.17/bwa";
our $SAMTOOLS = "/gscmnt/gc2737/ding/hsun/software/miniconda2/bin/samtools";
our $DISAMBIGUATE = "/gscmnt/gc2737/ding/hsun/software/miniconda2/bin/ngs_disambiguate";

our $HG_GENOME = "/gscmnt/gc2737/ding/hsun/data/human_genome/GRCh37-lite/GRCh37-lite.fa";
our $MM_GENOME = "/gscmnt/gc2737/ding/hsun/data/ensemble_v91/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa";


##-- Main
my @list_dir = `cat $folderList`;  #folder list

my $path = cwd;
`mkdir -p $path/job_log`;
my $record_log;    # record submitted job

my $name;
my $outScript;
foreach my $dir (@list_dir){
  chomp($dir);
  next if ($dir eq '');
	
  $dir =~ s/\/$//;
  $name = $1 if ($dir =~ /\/([.-\w]+)$/);
  
  $outScript = '';
  
  `mkdir -p $dir/log`;
  `mkdir -p $dir/script`;

	## bam2fq
  if ($program eq 'bam2fq'){
    $outScript = &Program_bam2fastq($dir, $name, $program);
	}
	  
  ## DNA-FQ-SortNameBam
  if ($program eq 'fq-sortBam' && $algorithm eq 'bwa'){
    $outScript = &Program_DNA_fqSortBam($dir, $name, $species, $program);
	}
	
	## Disambiguate
  if ($program eq 'disam' && $algorithm eq 'bwa'){
    $outScript = &Program_Disambiguate($dir, $name, $species, $program, $algorithm);
  }

	## New bam after Disambiguate
  if ($program eq 'newBam'){
    $outScript = &Program_afterDisambiguate($dir, $name, $program);
  }

	## DNA - Full step fq to Disambiguate
  if ($program eq 'fq2disam' && $algorithm eq 'bwa'){
    $outScript = &Program_DNA_fq2Disambiguate($dir, $name, $program, $algorithm);
  }  

	## RNA - Full step bam - fq - Disambiguate
  if ($program eq 'fq2disam' && $algorithm eq 'star'){
    $outScript = &Program_RNA_fq2Disambiguate($dir, $name, $program, $algorithm);
  }


  ## Run
	if (defined $run){
		if ($outScript ne ''){
			$record_log = &bsub2research_hpc("$dir/log", $name, $program, $outScript);
			sleep 1;
		} else {
			print $red, "WARNING: No matched script regarding ";
			print $normal, "$name.$program !\n\n";
		}
		
		# submit job record
		open my ($OUT),">$path/job_log/$name.$program.log";
		print $OUT "$record_log\n";
		close $OUT;
	}

}

exit;

##########################################################
## bsub to research-hpc
sub bsub2research_hpc
{
	my ($dir, $name, $program, $bashFile) = @_;
	
	my $lsf_log = "$dir/$name.$program.log";
	my $lsf_err = "$dir/$name.$program.err";

	`rm -f $lsf_log $lsf_err` if (-e $lsf_log || -e $lsf_err);
	
	# STAR software RAM 64GB >32GB
	my $bsub_cmd = "bsub -o $lsf_log -e $lsf_err -q research-hpc -M 64000000 -R \"select[mem>64000] rusage[mem=64000] span[hosts=1]\" -n 1 -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\' sh $bashFile\n";

	system($bsub_cmd);
	
	return($bsub_cmd);
}


## BamToFastQ
sub Program_bam2fastq
{
  my ($dir, $name, $program) = @_;
  
  my $outScript = "$dir/script/$name.$program.sh";
  
  open my ($OUT), '>', $outScript;
  
  print $OUT <<EOF	
#!/bin/bash

samtools=$SAMTOOLS

# bam2fq
\$samtools sort -m 1G -@ 6 -o $dir/$name.sortbyname.bam -n $dir/$name.bam
\$samtools fastq $dir/$name.sortbyname.bam -1 $dir/$name\_1.fastq.gz -2 $dir/$name\_2.fastq.gz
EOF
  ;

  close $OUT;
  
  return($outScript);

}



## DNA - fqSortByNameBam 
sub Program_DNA_fqSortBam
{
  my ($dir, $name, $species, $program) = @_;
  
  my $ref_genome;
  $ref_genome = $HG_GENOME if ($species eq 'human');
  $ref_genome = $MM_GENOME if ($species eq 'mouse');
    
  my $outScript = "$dir/script/$name.$program.$species.sh";
  
  open my ($OUT), '>', $outScript;
  
  print $OUT <<EOF	
#!/bin/bash

bwa=$BWA
samtools=$SAMTOOLS

genome=$ref_genome

# bwa hg/mm
\$bwa mem -t 8 -M -R "\@RG\\tID:$name\\tSM:$name\\tPL:illumina\\tLB:$name.lib\\tPU:$name.unit" \$genome $dir/$name\_1.fastq.gz $dir/$name\_2.fastq.gz > $dir/$species.sam
\$samtools view -Sbh $dir/$species.sam > $dir/$species.bam

# sort bam by natural name
\$samtools sort -m 1G -@ 6 -o $dir/$species.sortN.bam -n $dir/$species.bam

rm -f $dir/$species.sam $dir/$species.bam

# flagstat
\$samtools flagstat $dir/$species.sortN.bam > $dir/$species.sortN.bam.flagstat

EOF
  ;
    
  close $OUT;
  
  return($outScript);
}



## Disambiguate 
sub Program_Disambiguate
{
  my ($dir, $name, $species, $program, $algorithm) = @_;
  
  my $outScript = "$dir/script/$name.$program.sh";
  
  open my ($OUT), '>', $outScript;
  
  print $OUT <<EOF	
#!/bin/bash

disambiguate=$DISAMBIGUATE
samtools=$SAMTOOLS
bwa=$BWA
java=$JAVA
picard=$PICARD

hg_genome=$HG_GENOME

# mouse-filter - Disambiguate
\$disambiguate -s $name -o $dir -a $algorithm $dir/human.sortN.bam $dir/mouse.sortN.bam

# re-create fq
\$samtools sort -m 1G -@ 6 -o $dir/$name.disam.sortbyname.bam -n $dir/$name.disambiguatedSpeciesA.bam
\$samtools fastq $dir/$name.disam.sortbyname.bam -1 $dir/$name.disam_1.fastq.gz -2 $dir/$name.disam_2.fastq.gz

# bwa hg
\$bwa mem -t 8 -M -R "\@RG\\tID:$name\\tSM:$name\\tPL:illumina\\tLB:$name.lib\\tPU:$name.unit" \$hg_genome $dir/$name.disam_1.fastq.gz $dir/$name.disam_2.fastq.gz > $dir/$name.disam.reAligned.sam

# sort
\$java -Xmx8G -jar \$picard SortSam \\
   I=$dir/$name.disam.reAligned.sam \\
   O=$dir/$name.disam.reAligned.bam \\
   SORT_ORDER=coordinate

# remove-duplication
\$java -Xmx8G -jar \$picard MarkDuplicates \\
   I=$dir/$name.disam.reAligned.bam \\
   O=$dir/$name.disam.reAligned.remDup.bam \\
   REMOVE_DUPLICATES=true \\
   M=$dir/picard.$name.disam.reAligned.remdup.metrics.txt

# index bam
\$samtools index $dir/$name.disam.reAligned.remDup.bam

EOF
	;
  close $OUT;
  
  return($outScript);
}



## New bam after Disambiguate 
sub Program_afterDisambiguate
{
  my ($dir, $name, $program) = @_;
  
  my $outScript = "$dir/script/$name.$program.sh";
  
  open my ($OUT), '>', $outScript;
  
  print $OUT <<EOF	
#!/bin/bash

samtools=$SAMTOOLS
bwa=$BWA
java=$JAVA
picard=$PICARD

hg_genome=$HG_GENOME

# re-create fq
\$samtools sort -m 1G -@ 6 -o $dir/$name.disam.sortbyname.bam -n $dir/$name.disambiguatedSpeciesA.bam
\$samtools fastq $dir/$name.disam.sortbyname.bam -1 $dir/$name.disam_1.fastq.gz -2 $dir/$name.disam_2.fastq.gz

# bwa hg
\$bwa mem -t 8 -M -R "\@RG\\tID:$name\\tSM:$name\\tPL:illumina\\tLB:$name.lib\\tPU:$name.unit" \$hg_genome $dir/$name.disam_1.fastq.gz $dir/$name.disam_2.fastq.gz > $dir/$name.disam.reAligned.sam

# sort
\$java -Xmx8G -jar \$picard SortSam \\
   I=$dir/$name.disam.reAligned.sam \\
   O=$dir/$name.disam.reAligned.bam \\
   SORT_ORDER=coordinate

# remove-duplication
\$java -Xmx8G -jar \$picard MarkDuplicates \\
   I=$dir/$name.disam.reAligned.bam \\
   O=$dir/$name.disam.reAligned.remDup.bam \\
   REMOVE_DUPLICATES=true \\
   M=$dir/picard.$name.disam.reAligned.remdup.metrics.txt

# index bam
\$samtools index $dir/$name.disam.reAligned.remDup.bam

EOF
	;
  close $OUT;
  
  return($outScript);
}



## DNA-Full
sub Program_DNA_fq2Disambiguate
{
  my ($dir, $name, $program, $algorithm) = @_;
  
  my $outScript = "$dir/script/$name.$program.sh";
  
  open my ($OUT), '>', $outScript;
  
  print $OUT <<EOF	
#!/bin/bash

bwa=$BWA
samtools=$SAMTOOLS
disambiguate=$DISAMBIGUATE
java=$JAVA
picard=$PICARD

hg_genome=$HG_GENOME
mm_genome=$MM_GENOME

# bwa human
\$bwa mem -t 8 -M -R "\@RG\\tID:$name\\tSM:$name\\tPL:illumina\\tLB:$name.lib\\tPU:$name.unit" \$hg_genome $dir/$name\_1.fastq.gz $dir/$name\_2.fastq.gz > $dir/human.sam
\$samtools view -Sbh $dir/human.sam > $dir/human.bam

# sort bam by natural name
\$samtools sort -m 1G -@ 6 -o $dir/human.sortN.bam -n $dir/human.bam

rm -f $dir/human.sam $dir/human.bam

# bwa mouse
\$bwa mem -t 8 -M -R "\@RG\\tID:$name\\tSM:$name\\tPL:illumina\\tLB:$name.lib\\tPU:$name.unit" \$mm_genome $dir/$name\_1.fastq.gz $dir/$name\_2.fastq.gz > $dir/mouse.sam
\$samtools view -Sbh $dir/mouse.sam > $dir/mouse.bam

# sort bam by natural name
\$samtools sort -m 1G -@ 6 -o $dir/mouse.sortN.bam -n $dir/mouse.bam

rm -f $dir/mouse.sam $dir/mouse.bam

# mouse-filter - Disambiguate
\$disambiguate -s $name -o $dir -a $algorithm $dir/human.sortN.bam $dir/mouse.sortN.bam

# re-create fq
\$samtools sort -m 1G -@ 6 -o $dir/$name.disam.sortbyname.bam -n $dir/$name.disambiguatedSpeciesA.bam
\$samtools fastq $dir/$name.disam.sortbyname.bam -1 $dir/$name.disam_1.fastq.gz -2 $dir/$name.disam_2.fastq.gz

# bwa hg
\$bwa mem -t 8 -M -R "\@RG\\tID:$name\\tSM:$name\\tPL:illumina\\tLB:$name.lib\\tPU:$name.unit" \$hg_genome $dir/$name.disam_1.fastq.gz $dir/$name.disam_2.fastq.gz > $dir/$name.disam.reAligned.sam

# sort
\$java -Xmx8G -jar \$picard SortSam \\
   I=$dir/$name.disam.reAligned.sam \\
   O=$dir/$name.disam.reAligned.bam \\
   SORT_ORDER=coordinate

# remove-duplication
\$java -Xmx8G -jar \$picard MarkDuplicates \\
   I=$dir/$name.disam.reAligned.bam \\
   O=$dir/$name.disam.reAligned.remDup.bam \\
   REMOVE_DUPLICATES=true \\
   M=$dir/picard.$name.disam.reAligned.remdup.metrics.txt

# index bam
\$samtools index $dir/$name.disam.reAligned.remDup.bam

EOF
	;
  close $OUT;
  
  return($outScript);
}


## RNA-Full bam-fq-disambiguate
sub Program_RNA_fq2Disambiguate
{
  my ($dir, $name, $program, $algorithm) = @_;
  
  my $bam = "$dir/$name.bam"; # it useful when use bam2fq
  my $outFolder = $dir;       # for matching testing script

  my $outScript = "$dir/script/$name.$program.sh";
  open my ($OUT), '>', $outScript;
  
  print $OUT <<EOF
#!/bin/bash

star=/gscmnt/gc2737/ding/hsun/software/miniconda2/bin/STAR
samtools=/gscmnt/gc2737/ding/hsun/software/miniconda2/bin/samtools
disambiguate=/gscmnt/gc2737/ding/hsun/software/miniconda2/bin/ngs_disambiguate

hg_genomeDir=/gscmnt/gc2737/ding/hsun/data/ensemble_v91/GRCh37_star_genomeDir
hg_gtf=/gscmnt/gc2737/ding/hsun/data/ensemble_v91/Homo_sapiens.GRCh37.87.gtf
hg_genome=/gscmnt/gc2737/ding/hsun/data/ensemble_v91/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa

mm_genomeDir=/gscmnt/gc2737/ding/hsun/data/ensemble_v91/mm10_star_genomeDir
mm_gtf=/gscmnt/gc2737/ding/hsun/data/ensemble_v91/Mus_musculus.GRCm38.91.gtf
mm_genome=/gscmnt/gc2737/ding/hsun/data/ensemble_v91/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa


## bam2fq
## \$samtools sort -m 1G -@ 6 -o $outFolder/$name.sortbyname.bam -n $bam
## \$samtools fastq $outFolder/$name.sortbyname.bam -1 $outFolder/$name\_1.fastq.gz -2 $outFolder/$name\_2.fastq.gz

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
  close $OUT;
  
  return($outScript);
}





