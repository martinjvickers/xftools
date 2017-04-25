# methylation_tools

Getting this to compile on work machine;

```
~/cmake-3.6.0-rc2-Linux-x86_64/bin/cmake ../single-c_combine -DCMAKE_MODULE_PATH=../../seqan/util/cmake -DSEQAN_INCLUDE_PATH=../../seqan/include -DCMAKE_CXX_FLAGS=-std=c++14 -DCMAKE_BUILD_TYPE=Release
```

## Usage

Specify any number of input files using the `-i` flag.

```
single-c_combine -i file1.gff -i file2.gff ..... -i fileN.gff -o output.gff
```

To make a test file for diff

```
cat combined.single-c.gff | awk '{printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t%4f\t"$7"\t"$8"\t"$9"\n", $6}'
```
## Performance

### Original 

Can only take in two file at a time

```
gff_arithmetics.pl 'methyl()' -a -s d1d2_meiocyte_rep1.CG.tair10_sorted.gff d1d2_meiocyte_rep2.CG.tair10_sorted.gff > combined.single-c.gff
```

```
mvickers@x250 ~$ /usr/bin/time -v gff_arithmetics.pl 'methyl()' -a -s d1d2_meiocyte_rep1.CG.tair10_sorted.gff d1d2_meiocyte_rep2.CG.tair10_sorted.gff > combined.single-c.gff
	Command being timed: "gff_arithmetics.pl methyl() -a -s d1d2_meiocyte_rep1.CG.tair10_sorted.gff d1d2_meiocyte_rep2.CG.tair10_sorted.gff"
	User time (seconds): 1070.10
	System time (seconds): 18.34
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 18:11.49
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 40332
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 1
	Minor (reclaiming a frame) page faults: 231375
	Voluntary context switches: 14314
	Involuntary context switches: 60241
	Swaps: 0
	File system inputs: 873080
	File system outputs: 486480
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```

### Merging two files

```
mvickers@x250:~$ ls -lath example_data/d1d2_meiocyte_rep1.CHH.tair10.gff
-rwxr-xr-x 1 mvickers JIC_c1 1.3G Sep  6  2016 example_data/d1d2_meiocyte_rep1.CHH.tair10.gff
mvickers@x250:~$ ls -lath example_data/d1d2_meiocyte_rep2.CHH.tair10.gff
-rwxr-xr-x 1 mvickers JIC_c1 1.4G Sep  6  2016 example_data/d1d2_meiocyte_rep2.CHH.tair10.gff

mvickers@x250:~$ /usr/bin/time -v ./single-c_combine -i example_data/d1d2_meiocyte_rep1.CHH.tair10.gff -i example_data/d1d2_meiocyte_rep2.CHH.tair10.gff > meh
	User time (seconds): 220.30
	System time (seconds): 76.96
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 4:57.58
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 13019572
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 3254921
	Voluntary context switches: 61
	Involuntary context switches: 52907
	Swaps: 0
	File system inputs: 0
	File system outputs: 3281504
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```

### Merging three files

```
mvickers@x250:~$ /usr/bin/time -v ./single-c_combine -i example_data/d1d2_meiocyte_rep1.CHH.tair10.gff -i example_data/d1d2_meiocyte_rep2.CHH.tair10.gff -i example_data/d1d2.spm.CHH.tair10.gff > meh
ERROR: Unexpected end of input.
	User time (seconds): 371.37
	System time (seconds): 78.80
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 7:32.44
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 13096812
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 3274231
	Voluntary context switches: 86
	Involuntary context switches: 68118
	Swaps: 0
	File system inputs: 0
	File system outputs: 3350600
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```
