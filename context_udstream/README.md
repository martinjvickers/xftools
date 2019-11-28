## Context udstream


## Usage

```
context_udstream -r reference.fa -i sample_w1.gff -o out.gff
```

NOTE: The contig names need to be exactly the same (e.g. `Chr1` is not the same as `chr1`) in both the reference genome fasta file and the region GFF file.


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
    -s, --window-size INT
          Size of window Default: 50.
    -p, --percentage
          Rather than calculating the number of C's, calculate the percentage of C's in the window.
    -r, --input-region-file IN
          Path to the input file contains regions you're interested in calculating stats for.
```

## Complile

```
cmake ../context_udstream -DCMAKE_MODULE_PATH=../../seqan/util/cmake -DSEQAN_INCLUDE_PATH=../../seqan/include -DCMAKE_CXX_FLAGS=-std=c++14 -DCMAKE_BUILD_TYPE=Release
```


