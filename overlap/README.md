# methylation_tools

## Overlap Finder

This program is designed to find overlaps between a single-c w1.gff file and an annotation GFF file.

### Compiling

Getting this to compile on work machine;

```
~/cmake-3.6.0-rc2-Linux-x86_64/bin/cmake ../overlap -DCMAKE_MODULE_PATH=../../seqan/util/cmake -DSEQAN_INCLUDE_PATH=../../seqan/include -DCMAKE_CXX_FLAGS=-std=c++14 -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_ARGS="-DSEQAN_DISABLE_VERSION_CHECK=YES"
```

### Running

```
Overlap - Overlap
=================

SYNOPSIS
    Overlap -i input.gff -a annotation.gff -o output.gff [OPTIONS]

DESCRIPTION
    Given an annotation/feature file and an input w1 file, give counts for each feature

    -h, --help
          Display the help message.
    -i, --input-file IN
          Path to the input file
    -a, --input-annotation-file IN
          Path to the input filter file
    -o, --output-file OUT
          Path to the output file
    -l, --lazy-ref
          Internally it will capitalise both the input and annoation reference names so that chr1, Chr1 and CHR1 will
          all match. The output GFF will be of the same format as the annoation file.
    --version
          Display version information.
```

