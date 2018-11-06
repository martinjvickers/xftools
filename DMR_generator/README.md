# XFTOOLS

## DMR_generator

The purpose of this tool is to create "control" features of interest that mimick specific properties of "real" features however are randomly located across the genome.

The most basic of this is to create "control" features that have the same size distribution as the "real" features, as well as properties such as must overlap with gene, must be `X` distance from TE etc.

### Compiling

Getting this to compile on work machine;

```
~/cmake-3.6.0-rc2-Linux-x86_64/bin/cmake ../DMR_generator \
   -DCMAKE_MODULE_PATH=../seqan/util/cmake \
   -DSEQAN_INCLUDE_PATH=../seqan/include \
   -DCMAKE_CXX_FLAGS=-std=c++14 \
   -DCMAKE_BUILD_TYPE=Release \
   -DCMAKE_CXX_ARGS="-DSEQAN_DISABLE_VERSION_CHECK=YES" 
```

### Plan

I need the TE and Gene files as input.

I need the feature file to mimick.

I think it would be BEST to require the fasta file too.
   And then, if there is an issue where say the Gene annotation is CHR1-CHR5 but the fasta is CHR1-CHR5+CHRM+CHRC then the program will only return features from chromosomes found within the feature file. 

### Running

```
./DMR_generator -a genes.gff -f dmrs.gff -g genome.fa -l -e 200 -o mimicked.dmrs.gff
```


