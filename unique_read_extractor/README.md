# methylation_tools

## unique_read_extractor

This program extracts only unique reads from a tophat accepted hits file.

### Compiling

Getting this to compile on work machine;

```
~/cmake-3.6.0-rc2-Linux-x86_64/bin/cmake ../unique_read_extractor -DCMAKE_MODULE_PATH=../../seqan/util/cmake -DSEQAN_INCLUDE_PATH=../../seqan/include -DCMAKE_CXX_FLAGS=-std=c++14 -DCMAKE_BUILD_TYPE=Release
```

## Getting started

The program is run like this;

```
unique_read_extractor -i accepted_hits.sam -o unique.sam
```

You need to convert your BAM files into SAM files before hand;

```
samtools view -h -b -o accepted_hits.sam accepted_hits.bam
```

and then to convert the unique SAM file to a sorted BAM file for IGV;

```
samtools view -bS unique.sam | samtools sort - -o unique_sorted.bam
samtools index unique_sorted.bam
```

### MY ASSUMPTIONS SO FAR

* Since it's the tophat accepted hits file there are no unmapped reads
* Data is not sorted by read name hence the need to store names in RAM
* A unique read is defined as a read name that has only ever mapped once
