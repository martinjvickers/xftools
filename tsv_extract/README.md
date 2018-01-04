# XFTOOLS

## TOOLNAME

Description of TOOLNAME

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

### TODEL: Template instructions

TODEL: You should be able to `sed` for TOOLNAME within this doc as well as
TODEL: CMakeList.txt and TOOLNAME.cpp in order to rename. Also grep out
TODEL: the TODEL lines within this doc.
