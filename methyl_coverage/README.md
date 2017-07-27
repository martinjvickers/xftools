## How to do this

I need the three w1 files so that I can do the counting. I also need to know the number of C+G's within the reference genome.

./methyl_coverage --cg sample.cg.w1.gff --chg sample.chg.w1.gff --chh sample.chh.w1.gff -r ref.fasta


## Compile

```
~/cmake-3.6.0-rc2-Linux-x86_64/bin/cmake ../methyl_coverage -DCMAKE_MODULE_PATH=../../seqan/util/cmake -DSEQAN_INCLUDE_PATH=../../seqan/include -DCMAKE_CXX_FLAGS=-std=c++14 -DCMAKE_BUILD_TYPE=Release
```
