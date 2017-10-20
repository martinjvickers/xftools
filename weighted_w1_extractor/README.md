## weighted_w1_extractor

The purpose of this program is to create a GFF file containing the (weighted) number of times a base is covered by a read. A GFF which contains information about each base in our lab is known as a w1 file (AKA window of size 1). A regular w1 file is would simply count the number of times a base is covered by a read, if this is the functionality you desire, then use the bam_2_w1_extractor program. The weighted part of this is that if your alignment file contains reads which have mapped multiple times, e.g. for when looking at sRNA using bowtie2 searching and reporting all results, you may have a read that has mapped in several locations. If you have a read that maps to three locations, rather than counting that as +1 to each location, it will be counted as +(1/3) to each location. 

## Usage

```
weighted_w1_extractor -i test.bam -o test_w1.gff -l test_sample
```

## Options

```
    -h, --help
          Display the help message.
    --version
          Display version information.
    -i, --input-file IN
          Path to the input file
    -o, --output-file OUT
          Path to the output file
    -l, --label TEXT
          Column 3 GFF output label. Useful if using SignalMap as GFFs with the same label will be merged. Default:
          window.
```

## Complile

```
~/cmake-3.6.0-rc2-Linux-x86_64/bin/cmake ../weighted_w1_extractor -DCMAKE_MODULE_PATH=../seqan/util/cmake -DSEQAN_INCLUDE_PATH=../seqan/include -DCMAKE_CXX_FLAGS=-std=c++14 -DCMAKE_BUILD_TYPE=Release
```
