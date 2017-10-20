## Testing

Creating a combined file using the two input files

```
gff_arithmetics.pl 'methyl()' -a -s example_data/example_rep1.w1.gff example_data/example_rep2.w1.gff > old_program_result.w1_rep1_rep2.gff
```
The interesting part here is that there if you run this the other way around;

```
gff_arithmetics.pl 'methyl()' -a -s example_data/example_rep2.w1.gff example_data/example_rep1.w1.gff > old_program_result.w1_rep2_rep1.gff
```

you get different results

```
wc -l old_program_result.w1_rep1_rep2.gff old_program_result.w1_rep2_rep1.gff
  74615 old_program_result.w1_rep1_rep2.gff
  74622 old_program_result.w1_rep2_rep1.gff
 149237 total

```

So there is a bug somewhere in the original program. When rep1 is followed by rep2 in the input we appear to be missing several lines which appear in `example_rep1.w1.gff` but not in `example_rep2.w1.gff`;

```
diff old_program_result.w1_rep1_rep2.gff old_program_result.w1_rep2_rep1.gff
74615a74616,74622
> CHRM	.	CG	135939	135939	0.01	+	.	c=4;t=722
> CHRM	.	CG	135940	135940	0.00	-	.	c=1;t=563
> CHRM	.	CG	135944	135944	0.01	+	.	c=4;t=693
> CHRM	.	CG	135945	135945	0.00	-	.	c=1;t=587
> CHRM	.	CG	135954	135954	0.01	+	.	c=8;t=670
> CHRM	.	CG	135955	135955	0.00	-	.	c=0;t=590
> CHRM	.	CG	135966	135966	0.01	+	.	c=9;t=682
```

To ensure this is correct, I first check to see if the results are the same both ways around and then I compare the results file against the working way around (`example_data/old_program_result.w1_rep2_rep1.gff`).

### A note on diff

The perl and C++ code differ on the way they print the score data. So in order to simply run `diff` between my output and the output from the old program I have used printf to force the score to be a float at 2 decimal places.

```
cat result.gff | awk '{printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t%.2f\t"$7"\t"$8"\t"$9"\n", $6}' > result_mod.gff
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
