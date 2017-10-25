## cytosine_abundance

The purpose of this program is to calculate the cytosine abundance and cytosine context's in a reference genome.

## Usage

```
cytosine_abundance -i reference.fa -o reference_cytosine.gff
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
~/cmake-3.6.0-rc2-Linux-x86_64/bin/cmake ../cytosine_abundance -DCMAKE_MODULE_PATH=../seqan/util/cmake -DSEQAN_INCLUDE_PATH=../seqan/include -DCMAKE_CXX_FLAGS=-std=c++14 -DCMAKE_BUILD_TYPE=Release
```
