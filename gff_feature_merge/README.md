# methylation_tools

## GFF Feature Merge

### Compiling

Getting this to compile on work machine;

```
~/cmake-3.6.0-rc2-Linux-x86_64/bin/cmake ../gff_feature_merge -DCMAKE_MODULE_PATH=../../seqan/util/cmake -DSEQAN_INCLUDE_PATH=../../seqan/include -DCMAKE_CXX_FLAGS=-std=c++14 -DCMAKE_BUILD_TYPE=Release
```

### Getting started

```
./gff_feature_merge -i example_data/reads_only_in_iae.gff -s 2 -o meh.txt 
```

###What does it do?

This program simply merges elements/features of a GFF file that are close to oneanother given a user defined variable (`-s`). The best way to explain this is by and example. Consider the following input GFF file.

```
chr2	.	window	80	80	1	.	.	n=1
chr2	.	window	81	81	1	.	.	n=1
chr2	.	window	82	82	1	.	.	n=1
chr2	.	window	90	90	1	.	.	n=1
chr2	.	window	91	91	1	.	.	n=1
chr2	.	window	101	101	1	.	.	n=1
chr2	.	window	102	102	1	.	.	n=1
chr2	.	window	104	104	1	.	.	n=1
chr2	.	window	200	200	1	.	.	n=1
chr2	.	window	600	600	1	.	.	n=1
```

The aim of this program is to merge features (lines/rows) within the GFF that are within `-s` of oneanother. So, if `-s 2` is applied, then the resulting output will be;

```
chr2    .       window  80      82      .       .       .       n=1
chr2    .       window  90      91      .       .       .       n=1
chr2    .       window  101     104     .       .       .       n=1
chr2    .       window  200     200     .       .       .       n=1
chr2    .       window  600     600     .       .       .       n=1
```

similarly if `-s 10` is applied the resulting output will be;

```
chr2    .       window  80      104      .       .       .       n=1
chr2    .       window  200     200     .       .       .       n=1
chr2    .       window  600     600     .       .       .       n=1
```

