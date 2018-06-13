# XFTOOLS

## max_wba

Given an annotation/feature file and an input GFF file, it will return the record in the input GFF with the
maximum score that resides within the range of the annotation file. If more than one record is found with equal
max score, it returns each with that max score. NOTE: this program is new and relatively untested.

### Compiling

Getting this to compile on work machine;

```
~/cmake-3.6.0-rc2-Linux-x86_64/bin/cmake ../max_wba \
   -DCMAKE_MODULE_PATH=../../seqan/util/cmake \
   -DSEQAN_INCLUDE_PATH=../../seqan/include \
   -DCMAKE_CXX_FLAGS=-std=c++14 \
   -DCMAKE_BUILD_TYPE=Release \
   -DCMAKE_CXX_ARGS="-DSEQAN_DISABLE_VERSION_CHECK=YES" 
```

### Running

```
$ max_wba -i example_data/input_file.w50.gff -a example_data/annotation_file.dmrs.gff -o example_data/results.w50.gff
```


```
$ max_wba --help
max_wba - WBA, returning max score record
=========================================

SYNOPSIS
    max_wba -i input.w50.gff -a annotation.dmrs.gff -o output.gff [OPTIONS]

DESCRIPTION
    Given an annotation/feature file and an input GFF file, it will return the record in the input GFF with the
    maximum score that resides within the range of the annotation file. If more than one record is found with equal
    max score, it returns each with that max score. NOTE: this program is new and relatively untested

OPTIONS
    -h, --help
          Display the help message.
    --version-check BOOL
          Turn this option off to disable version update notifications of the application. One of 1, ON, TRUE, T, YES,
          0, OFF, FALSE, F, and NO. Default: 1.
    -i, --input-file INPUT_FILE
          Input data file, e.g. w50s
    -a, --input-annotation-file INPUT_FILE
          Input annotation file, e.g. DMRs
    -o, --output-file OUTPUT_FILE
          Path to the output file
    --version
          Display version information.

VERSION
    Last update: June 2018
    max_wba version: 0.0.1
    SeqAn version: 2.3.2
```

### Detailed purpose

The idea, given a DMR/annotation file of some sort along with a second data file (e.g. w50 or w1), return the result with
the maximum score (column 6) which overlaps the annotation.

Annotation/DMR file:

```
CHR1    .       WT.SN.SHXF18A-WT.VN.SHXF2H.tair10.w50.CG.dif0.5p0.001.merged.d200bp.500_1500bp  193350  193900  .       .       .       ID=DMR_1
CHR1    .       WT.SN.SHXF18A-WT.VN.SHXF2H.tair10.w50.CG.dif0.5p0.001.merged.d200bp.500_1500bp  422750  423300  .       .       .       ID=DMR_2
CHR1    .       WT.SN.SHXF18A-WT.VN.SHXF2H.tair10.w50.CG.dif0.5p0.001.merged.d200bp.500_1500bp  432300  432950  .       .       .       ID=DMR_3
CHR1    .       WT.SN.SHXF18A-WT.VN.SHXF2H.tair10.w50.CG.dif0.5p0.001.merged.d200bp.500_1500bp  639700  641200  .       .       .       ID=DMR_4
```

Data file (in this case w50 file):

```
CHR1    .       WT.SN.SHXF18A-WT.VN.SHXF2H.tair10.w50.CG.dif0.5p0.001   193351  193400  .55309999       .       .       p=5.3510398e-14
CHR1    .       WT.SN.SHXF18A-WT.VN.SHXF2H.tair10.w50.CG.dif0.5p0.001   193551  193600  .57599998       .       .       p=2.3102299e-07
CHR1    .       WT.SN.SHXF18A-WT.VN.SHXF2H.tair10.w50.CG.dif0.5p0.001   193601  193650  .69489998       .       .       p=7.8985698e-18
CHR1    .       WT.SN.SHXF18A-WT.VN.SHXF2H.tair10.w50.CG.dif0.5p0.001   193651  193700  .9285   .       .       p=0
CHR1    .       WT.SN.SHXF18A-WT.VN.SHXF2H.tair10.w50.CG.dif0.5p0.001   193851  193900  .61080003       .       .       p=1.5083300e-10
```

Result

```
CHR1    .       WT.SN.SHXF18A-WT.VN.SHXF2H.tair10.w50.CG.dif0.5p0.001   193651  193700  0.9285  .       .       p=0
CHR1    .       WT.SN.SHXF18A-WT.VN.SHXF2H.tair10.w50.CG.dif0.5p0.001   423151  423200  0.8948  .       .       p=5.2231999e-15
CHR1    .       WT.SN.SHXF18A-WT.VN.SHXF2H.tair10.w50.CG.dif0.5p0.001   432601  432650  0.861   .       .       p=1.7892700e-29
CHR1    .       WT.SN.SHXF18A-WT.VN.SHXF2H.tair10.w50.CG.dif0.5p0.001   640551  640600  0.8382  .       .       p=3.0714699e-28
```

NOTE: If two or more entries in the data/w50 file have the same score, all of those with the same max score will be returned.

### Notes

* This was done super quick for Shengbo to get his work done, so it's not tested
* Uses a brute force loop-of-loop method to finish this so VERY slow
* I need to rewrite when I have more time, implementing my WBA algorithm.
