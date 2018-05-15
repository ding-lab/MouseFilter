## MouseFilter

### GMT Filter-mouse-bases
```
NOTE:Please refer to example.mgi.gmt.sh
     If test the example.mgi.gmt.sh, you must modify the vars including name, outDir and bed.
     The shell script only run in MGI servers.
     
     The *.permissive.out is the final filtered result.
```

### Disambiguate

* For running mouse filtering in research-hpc, please use "createBash.disambiguate.v2.pl".

* The "run.lsf.disambiguate.v2.1.pl" is the testing version.
The script can filter mouse reads of DNA- and RNA-seq NGS.
Created the script purpose is to test and create CWL script.

     NOTE: The MGI -q long is a unstable queue now.

* Dockers (https://github.com/ding-lab/dockers)
    * Only Disambiguate
   
          docker pull hsun9/disambiguate
          docker run hsun9/disambiguate ngs_disambiguate --help
    
    * Full pipeline of WXS docker image
    
          docker pull hsun9/disambiguateplus


