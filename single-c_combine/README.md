# single-c_combine

Merge several single-c (AKA w1/window 1 files) together. This can be useful when you have had sequencing done, mapped and extracted the information you want (e.g. w1/w50) and then had more sequencing done to boost coverage.

## Usage

Specify any number of input files using the `-i` flag.

```
single-c_combine -i file1.gff -i file2.gff ..... -i fileN.gff -o output.gff
```

## Complile

Getting this to compile on work machine;

```
~/cmake-3.6.0-rc2-Linux-x86_64/bin/cmake ../single-c_combine -DCMAKE_MODULE_PATH=../../seqan/util/cmake -DSEQAN_INCLUDE_PATH=../../seqan/include -DCMAKE_CXX_FLAGS=-std=c++14 -DCMAKE_BUILD_TYPE=Release
```
