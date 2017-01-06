# methylation_tools

## FAST Extract

A simple tool to extract reads from a fasta/fastq (zipped or unzipped) given a list of IDs in the form of a text file. 

This could be done with GREP, but the key feature is that you can also exclude reads by ID from your input fast{aq} file.

### Compiling

Getting this to compile on work machine;

```
~/cmake-3.6.0-rc2-Linux-x86_64/bin/cmake ../fast_extract -DCMAKE_MODULE_PATH=../../seqan/util/cmake -DSEQAN_INCLUDE_PATH=../../seqan/include -DCMAKE_CXX_FLAGS=-std=c++14 -DCMAKE_BUILD_TYPE=Release
```

### Getting started

Simply download the static binary release and run the program. There are two modes to this.

To produce a file containing only reads whose IDs exist in the list_to_extract.txt file.

```
./fast_extract -i rawreads.fastq -f list_to_extract.txt > extracted_reads.fastq
```

If you wish to do the opposite, exlude (-e) the reads whose IDs exist in the list_to_extract.txt file, include the -e flag.

```
./fast_extract -e -i rawreads.fastq -f list_to_extract.txt > extracted_reads.fastq
```

### Preparing a list.txt file of IDs

This is the tricky bit I guess, but should be relatively simple if you understand how fast_extract matches IDs using SeqAn. In our example data a raw read looks like;

```
@K00114:151:H3C3JBBXX:2:1101:7140:1773 1:N:0:AAGCCTAT
GGGCAAGAGCCAGGCCTCGATGAGTAGGAGGGCGCGGCGGTCGCTGCA
+
<<FFFKA<FAK<AFF(AA7AKK<<KKAAAFF<A<F<FKFK<7AF,AA7
```

and entries in the list look like;

```
K00114:151:H3C3JBBXX:2:1222:28151:5219
K00114:151:H3C3JBBXX:2:2201:7293:36670
K00114:151:H3C3JBBXX:2:2128:13677:35035
K00114:151:H3C3JBBXX:2:2128:13535:35282
K00114:151:H3C3JBBXX:2:2208:24659:31168
K00114:151:H3C3JBBXX:2:2223:5567:45250
```

There are two things you need to know when preparing your list file. 
* Do not include the @ symbol at the beginning of the ID in your list file.
* Only include up until the first space or tab. 

So if you wish to extract the read above, the ID in the list file should be;

```
K00114:151:H3C3JBBXX:2:1101:7140:1773
```

### Performance

## Using a Map

```
mvickers@n108379:~/development/methylation_tools/fast_extract$ time `./fast_extract -e -i example_data/rawreads.fastq -f example_data/list_to_extract.txt > meh.fq`
real	0m0.868s
real	0m0.893s
real	0m0.867s
```
```
mvickers@n108379:~/development/methylation_tools/fast_extract$ time `./fast_extract -i example_data/rawreads.fastq -f example_data/list_to_extract.txt > meh.fq`
real	0m0.344s
real	0m0.346s
real	0m0.349s
```

## Using an Unordered Map

```
mvickers@n108379:~/development/methylation_tools/fast_extract$ time `./fast_extract -e -i example_data/rawreads.fastq -f example_data/list_to_extract.txt > meh.fq`
real    0m0.680s
real    0m0.683s
real    0m0.675s
```

```
mvickers@n108379:~/development/methylation_tools/fast_extract$ time `./fast_extract -i example_data/rawreads.fastq -f example_data/list_to_extract.txt > meh.fq`
real	0m0.144s
real	0m0.148s
real	0m0.156s
```
