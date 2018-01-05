# XFTOOLS

## tsv_extract

This tool is specifically written to process the Arabidopsis 1001 data. It's not a generic tool.
Do not use this on anything other than what it was designed for, specifically the TSV files from;

(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43857)

### Compiling

Getting this to compile on work machine;

```
~/cmake-3.6.0-rc2-Linux-x86_64/bin/cmake ../tsv_extract \
   -DCMAKE_MODULE_PATH=../../seqan/util/cmake \
   -DSEQAN_INCLUDE_PATH=../../seqan/include \
   -DCMAKE_CXX_FLAGS=-std=c++14 \
   -DCMAKE_BUILD_TYPE=Release \
   -DCMAKE_CXX_ARGS="-DSEQAN_DISABLE_VERSION_CHECK=YES" 
```

### Running

To run, pipe your tsv file into the program, while defining your annotation file and a sample prefix 
for the output files. This will give you four files, three gffs with methylation scores calculated
and a txt file with the raw information segmented out of the tsv.

```
$ zcat example_data/GSM1085191_mC_calls_Aa_0.tsv.gz | ./tsv_extract -a example_data/annotation_mod.gff -p samplename
```

