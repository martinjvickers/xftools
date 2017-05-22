#methylation_tools

Getting this to compile on work machine;

```
~/cmake-3.6.0-rc2-Linux-x86_64/bin/cmake ../w50_creator -DCMAKE_MODULE_PATH=../../seqan/util/cmake -DSEQAN_INCLUDE_PATH=../../seqan/include -DCMAKE_CXX_FLAGS=-std=c++14 -DCMAKE_BUILD_TYPE=Release
```

## Usage

Usage for methylation data is simple. Download the latest release and extract;

e.g.


```
wget https://github.com/martinjvickers/xftools/releases/download/v0.0.3/methylation_tools_v0.0.3.tar.gz .
tar xvf methylation_tools_v0.0.3.tar.gz
```

If you wish to try some example data;

```
wget https://github.com/martinjvickers/methylation_tools/blob/master/w50_creator/example_data/GSM952438_TFH39.all.cg-col.w50.gff.gz
```

and then here you can run the program which will bin the methyalation data from each single-c (w1) location in the example file into 50bp windows (-s flag). Since the -l flag is set the third column in the output is labeled as desired. This is useful for importing into SignalMap.

```
./w50_creator -s 50 -l GSM952438_CG -i GSM952438_TFH39.all.cg-col.gff.gz > GSM952438_TFH39.all.w50.gff
```

You can also use this for counts data by specifying the `-t count` flag which has different behaviour than the default `-t methyl` flag. `-t count` will sum columns 6 of rows that belong in a bin and also display a count in column 9, e.g. `n=10`. The command for this is as follows;

```
./w50_creator -i example_data/reads_start_bases.gff -t count -s 1 > reads_summed_w1.gff
```

## Testing

### Test Machine
Intel(R) Core(TM) i5-4590 CPU @ 3.30GHz
32GB RAM


### Input file:
```
-rw-r--r--  1 mvickers JIC_c1  32M Aug 15 13:44 GSM952438_TFH39.all.cg-col.gff
mvickers@n108379:~/development/w50_creator$ wc -l GSM952438_TFH39.all.cg-col.gff
663526 GSM952438_TFH39.all.cg-col.gff
```

### Perl script - old way
```
User time (seconds): 43.50
Maximum resident set size (kbytes): 112268
```

### My program - new way
```
User time (seconds): 1.13
Maximum resident set size (kbytes): 5980
```

### Input file:
```
-rw-------  1 mvickers JIC_c1 779M Aug 15 14:16 D73_rep3_PE_SRR521046_R1.fastq.gz_bismark_bt2_pe.CHH.w1.gff
mvickers@n108379:~/development/w50_creator$ wc -l D73_rep3_PE_SRR521046_R1.fastq.gz_bismark_bt2_pe.CHH.w1.gff
19743126 D73_rep3_PE_SRR521046_R1.fastq.gz_bismark_bt2_pe.CHH.w1.gff
```

### Perl script - old way
```
User time (seconds): 316.19
Maximum resident set size (kbytes): 12447076
```

### My program - new way
```
User time (seconds): 23.50
Maximum resident set size (kbytes): 194288
```

### Input file:
```
-rw-r--r--  1 mvickers JIC_c1  14G Jul 10 22:17 70SC_1_val_1.fq.gz_bismark_bt2_pe.CHH.w1.gff
mvickers@n108379:/mnt/ssd/mvickers/massive_speed_test$ wc -l 70SC_1_val_1.fq.gz_bismark_bt2_pe.CHH.w1.gff
358238974 70SC_1_val_1.fq.gz_bismark_bt2_pe.CHH.w1.gff
```

### My program - new way
```
User time (seconds): 415.12
Maximum resident set size (kbytes): 308248
```
