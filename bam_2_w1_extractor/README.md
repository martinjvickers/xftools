## Compile

Only seems to work with a later version of seqan. I don't know why.

```
~/cmake-3.6.0-rc2-Linux-x86_64/bin/cmake ../bam_2_w1_extractor -DCMAKE_MODULE_PATH=../seqan/util/cmake -DSEQAN_INCLUDE_PATH=../seqan/include -DCMAKE_CXX_FLAGS=-std=c++14 -DCMAKE_BUILD_TYPE=Release
```
