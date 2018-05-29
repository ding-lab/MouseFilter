
## Description
* __INPUT__ 
     
     *.bam or *.fq.gz (pair-end reads)
     
* __OUTPUT__

     WXS/WGS: output human *.bam & *.bai files, which filtered mouse reads, sorted and removed duplicated reads files.
     RNA-seq: output human *.fq.gz, which filtered mouse reads files.


## Install software
   * Java (JRE1.8.x)   # for running in MGI server
   * conda install -c bioconda bwa=0.7.15
   * conda install -c bioconda samtools=1.5
   * conda install -c bioconda picard=2.17.11
   * conda install -c bioconda bamtools=2.4.0
   * conda install -c bioconda ngs-disambiguate=2016.11.10
   * conda install -c bioconda star=2.5.4a

-------------------------------------------------------------------

## MouseFilter: Docker version

* Docker files (https://github.com/ding-lab/dockers)
* Usage
     
     ```
      * Only Disambiguate
   
          docker pull hsun9/disambiguate
          docker run hsun9/disambiguate ngs_disambiguate --help


      * Full pipeline of mouse filter (wxs data) docker image

          docker pull hsun9/disambiguateplus
          docker run hsun9/disambiguateplus ngs_disambiguate --help
     ```

-------------------------------------------------------------------

## MouseFilter: CWL version

The CWL version mainly developed by Matthew Wyczalkowski
(https://github.com/ding-lab/MouseTrap2)
     
-------------------------------------------------------------------

## MouseFilter: MGI-server version

NOTE : The pipeline tested in MGI cluster. 

### I) VCF/BED file based mouse filtering (gmt somatic filter-mouse-bases):
```
NOTE:Please refer to "example.mgi.gmt.sh"
     If test the example.mgi.gmt.sh, you must modify the vars including name, outDir and bed.
     The shell script only run in MGI servers.
     
     The *.permissive.out is the final filtered result.
```

```
name="example"                   #-- any name
outDir=/example/outFolder        #-- any dir
bed=/example/target.bed          #-- any bed file
hg19=/path/GRCh37-lite.fa        #-- GRCh37
mm10=/path/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa
chain10=/path/hg19ToMm10.over.chain

gmt somatic filter-mouse-bases --chain-file=$chain10 \
    --human-reference=$hg19 \
    --mouse-reference=$mm10 --variant-file=$bed \
    --filtered-file=$outDir/log/$name.mouse.hg19toMm10.out \
    --output-file=$outDir/$name.hg19toMm10.permissive.out --permissive

## target.bed (create based on VCF)
1 121009 121009 C T

## chain download 
http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToMm10.over.chain.gz
```



### II) FASTQ/BAM based mouse reads filtering (Disambiguate tool):
* For running mouse filtering in research-hpc, please use "createBash.disambiguate.v2.pl".

* The "run.lsf.disambiguate.v2.1.pl" is the testing version.
The script can filter mouse reads of DNA- and RNA-seq NGS.
Created the script purpose is to test and create CWL script.

     NOTE: The MGI -q long is a unstable queue now.


* __Mouse reads filter from DNA-seq based data (WGS/WXS)__
    ```
    # Func. Filter mouse reads and create sorted bam for running somatic calling using pair-end fastq
    # check created bash file
    perl createBash.disambiguate.v2.pl -f folderList -p fq2disam -a bwa -t dna
    
    # submit job to research-hpc
    perl createBash.disambiguate.v2.pl -f folderList -p fq2disam -a bwa -t dna -run
     ```
* __Mouse reads filter from RNA-seq based data__
     ```
     # Func. Filter mouse reads and create new fastq using pair-end fastq 
     # check created bash file
     perl createBash.disambiguate.v2.pl -f folderList -p fq2disam -a star
     
     # submit job to research-hpc
     perl createBash.disambiguate.v2.pl -f folderList -p fq2disam -a star -t rna -run
     ```
    
    ```
     # Role of folder
      e.g.  
          /data/240_wxs/240_wxs_1.fastq.gz
          /data/240_wxs/240_wxs_1.fastq.gz
          /data/240_wxs/240_wxs.bam

     # FolderList
     e.g.
          /data/240_wxs
          /data/241_wxs

     # Download reference:
     GRCh37 (or GRCh38) and GRCm38 from Ensembl
     For RNA-seq filter, GTF file must download
    
     NOTE: Before running perl script, the reference and software directory should change.
     ```


     
     
