# methylation_tools

Getting this to compile on work machine;

```
~/cmake-3.6.0-rc2-Linux-x86_64/bin/cmake ../w50_creator -DCMAKE_MODULE_PATH=../../seqan/util/cmake -DSEQAN_INCLUDE_PATH=../../seqan/include -DCMAKE_CXX_FLAGS=-std=c++14 -DCMAKE_BUILD_TYPE=Release
```

##Testing

###Test Machine
Intel(R) Core(TM) i5-4590 CPU @ 3.30GHz
32GB RAM


###Input file:
```
-rw-r--r--  1 mvickers JIC_c1  32M Aug 15 13:44 GSM952438_TFH39.all.cg-col.gff
mvickers@n108379:~/development/w50_creator$ wc -l GSM952438_TFH39.all.cg-col.gff
663526 GSM952438_TFH39.all.cg-col.gff
```

###Perl script - old way
```
User time (seconds): 43.50
Maximum resident set size (kbytes): 112268
```

###My program - new way
```
User time (seconds): 1.13
Maximum resident set size (kbytes): 5980
```

###Input file:
```
-rw-------  1 mvickers JIC_c1 779M Aug 15 14:16 D73_rep3_PE_SRR521046_R1.fastq.gz_bismark_bt2_pe.CHH.w1.gff
mvickers@n108379:~/development/w50_creator$ wc -l D73_rep3_PE_SRR521046_R1.fastq.gz_bismark_bt2_pe.CHH.w1.gff
19743126 D73_rep3_PE_SRR521046_R1.fastq.gz_bismark_bt2_pe.CHH.w1.gff
```

###Perl script - old way
```
User time (seconds): 316.19
Maximum resident set size (kbytes): 12447076
```

###My program - new way
```
User time (seconds): 23.50
Maximum resident set size (kbytes): 194288
```

