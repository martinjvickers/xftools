# Pausing Index

## Purpose

A program to take the TAIR10 exons as a GFF (this looks like this);

```
Chr1	TAIR10	exon	3631	3913	.	+	.	Parent=AT1G01010.1
Chr1	TAIR10	exon	3996	4276	.	+	.	Parent=AT1G01010.1
Chr1	TAIR10	exon	4486	4605	.	+	.	Parent=AT1G01010.1
Chr1	TAIR10	exon	4706	5095	.	+	.	Parent=AT1G01010.1
Chr1	TAIR10	exon	5174	5326	.	+	.	Parent=AT1G01010.1
Chr1	TAIR10	exon	5439	5899	.	+	.	Parent=AT1G01010.1
Chr1	TAIR10	exon	8571	8737	.	-	.	Parent=AT1G01020.1
Chr1	TAIR10	exon	8417	8464	.	-	.	Parent=AT1G01020.1
Chr1	TAIR10	exon	8236	8325	.	-	.	Parent=AT1G01020.1
Chr1	TAIR10	exon	7942	7987	.	-	.	Parent=AT1G01020.1
```

And a `w1` file, and then create a GFF labeling the TSS region based upon the start of the gene and the TSS size parameter `--tss-size`. 

## Install

Getting this to compile on work machine;

```
~/cmake-3.6.0-rc2-Linux-x86_64/bin/cmake ../pausing_index -DCMAKE_MODULE_PATH=../../seqan/util/cmake -DSEQAN_INCLUDE_PATH=../../seqan/include -DCMAKE_CXX_FLAGS=-std=c++14 -DCMAKE_BUILD_TYPE=Release
```

## Usage

```
./pausing_index -i example_data/example.w1.gff -a example_data/example_annotation.gff -tss 200 -o meh.gff
```
