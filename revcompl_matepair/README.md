# XFTOOLS

## reverse_complement

A simple program for reversing the complement of a fastq file. This is especially useful if, for example, you wish to run paired end BS reads through a directional BS mapper as though they are single-end reads. The first mate pair (usually R1 on an illumina machine) would be left alone, however the second mate pair (R2) would need to be reverse complemented.

### Compiling

Getting this to compile on work machine;

```
~/cmake-3.6.0-rc2-Linux-x86_64/bin/cmake ../revcompl_matepair \
   -DCMAKE_MODULE_PATH=../seqan/util/cmake \
   -DSEQAN_INCLUDE_PATH=../seqan/include \
   -DCMAKE_CXX_FLAGS=-std=c++14 \
   -DCMAKE_BUILD_TYPE=Release \
   -DCMAKE_CXX_ARGS="-DSEQAN_DISABLE_VERSION_CHECK=YES" 
```

### Running

```
revcompl_matepair -i input.fq.gz -o output.fq.gz
```


