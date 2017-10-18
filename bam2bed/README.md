## What and why

This is a really simple program to convert reads into a bed file. This is VERY basic right now due to time constraints.

Why would I want this? I needed to cluster the sRNA reads. The quickest program I knew to do this is `bedtools merge`, so rather than writing something to do the whole lot, I just wanted to convert the BAM/SAM file into a BED file so that `bedtools merge` would understand it.

## What about the CIGAR score?

Good question, this is why I used SeqAn rather than just ```awk '{print $3"\t"$4"\t"$4+length($10)}'```. The `-c` flag in the example below takes into account the CIGAR score.

NOTE: This will produce a BED record that is from the start to the end of the read where it maps in the reference, so includes insertions/deletions. But it does not deal with gaps, so if you have to insert in the reference (e.g. an `I`). e.g.

```
POS	12345
REF	TT-GT
	|| ||
READ	TTAGT
```

So, in this case, the start position would be 1 and the end position would be 5 with respect to the reference. 

For our purposes this is what we want, but it may not be what you want so be careful.

## Usage

```
./bam2bed -i reads.bam -o meh.bed -c > output.bed
```

NOTE: There is a sacrifical -o file right now but for some odd reason SeqAn was writing out the BED file oddly. So for now I am just printing to screen. I will get around to fixing this.

## Compile

```
~/cmake-3.6.0-rc2-Linux-x86_64/bin/cmake ../bam2bed -DCMAKE_MODULE_PATH=../seqan/util/cmake -DSEQAN_INCLUDE_PATH=../seqan/include -DCMAKE_CXX_FLAGS=-std=c++14 -DCMAKE_BUILD_TYPE=Release
make -j4
```
