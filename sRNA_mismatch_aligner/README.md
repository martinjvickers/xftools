## sRNA_mismatch_aligner


## Usage

```
sRNA_mismatch_aligner -i reference.fa -o reference_cytosine.gff
```


## Complile

```
~/cmake-3.6.0-rc2-Linux-x86_64/bin/cmake ../sRNA_mismatch_aligner -DCMAKE_MODULE_PATH=../seqan/util/cmake -DSEQAN_INCLUDE_PATH=../seqan/include -DCMAKE_CXX_FLAGS=-std=c++14 -DCMAKE_BUILD_TYPE=Release
```

