# methylation_tools

## FAST Extract

A simple tool to extract reads from a fasta/fastq (zipped or unzipped) given a list of IDs in the form of a text file. 

This could be done with GREP, but the key feature is that you can also exclude reads by ID from your input fast{aq} file.

### Compiling

Getting this to compile on work machine;

```
~/cmake-3.6.0-rc2-Linux-x86_64/bin/cmake ../fast_extract -DCMAKE_MODULE_PATH=../../seqan/util/cmake -DSEQAN_INCLUDE_PATH=../../seqan/include -DCMAKE_CXX_FLAGS=-std=c++14 -DCMAKE_BUILD_TYPE=Release
```

### Getting started

Simply download the static binary release and run as follows;

```
./fast_extract -e -i rawreads.fastq -f list_to_extract.txt > extracted_reads.fastq
```
