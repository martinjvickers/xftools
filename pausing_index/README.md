# Pausing Index

## Purpose

Given a GFF containing single base counts (i.e. from RNA-Seq) and a GFF containing exons this program will calculate the ratio between the rpkm at TSS and the rpkm for the rest of the coding exons. Known as the "pausing index".

## Usage

Quick start:

```
./pausing_index -i example_data/example.w1.gff -a example_data/example_annotation.gff -tss 200 -o out.gff
```



```
pausing_index - XFTOOLS
=======================

SYNOPSIS
    pausing_index -i input.w1.gff -a reference.gff -tss 200 [OPTIONS]

DESCRIPTION
    Calculates the pausing index of genes

OPTIONS
    -h, --help
          Display the help message.
    --version-check BOOL
          Turn this option off to disable version update notifications of the application. One of 1, ON, TRUE, T, YES,
          0, OFF, FALSE, F, and NO. Default: 1.
    -i, --input-file INPUT_FILE
          Path to the input file
    -a, --annotation-file INPUT_FILE
          Path to the input file
    -o, --output-file INPUT_FILE
          Path to the input file
    --version
          Display version information.
    -tss, --tss-size INTEGER
          Size of the TSS you wish to extract Default: 200.
```

## Files

### Annotation

The annotation file `-a` should only contain the exons and therefore look like this;

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

### Input

The input file `-i` should contain our `w1` file, a GFF containing single bases

```
chr1	xftools	window	5978	5978	1	.	.	n=1
chr1	xftools	window	6000	6000	2	.	.	n=2
chr1	xftools	window	6017	6017	1	.	.	n=1
chr1	xftools	window	6020	6020	1	.	.	n=1
chr1	xftools	window	6063	6063	3	.	.	n=3
chr1	xftools	window	6064	6064	1	.	.	n=1
chr1	xftools	window	6075	6075	2	.	.	n=2
chr1	xftools	window	6192	6192	1	.	.	n=1
chr1	xftools	window	6193	6193	2	.	.	n=2
chr1	xftools	window	6194	6194	1	.	.	n=1
```

It is the score in column 6 that is used for the calculations.

### Input File Compatibility

You may notice above that the references are different, `Chr1` and `chr1`. This is fine as internally `pausing_index` captilises everything. You just need to ensure that in this case the reference is `CHR1`, `Chr1` or something similar. You couldn't mix `1` with `chr1`. 

### Output file

The `-o` flag specifies the output file. This will also be a GFF and will appear as so;

```
Chr1    pausing_index   gene    710019  711803  .       -       .       Parent=AT1G03055.1;SumExonLenth=919;SumTSS=0.000000;SumGene=0.000000;rpkmTSS=0.000000;rpkmGene=0.000000
Chr1    pausing_index   gene    710064  711803  .       -       .       Parent=AT1G03055.2;SumExonLenth=1027;SumTSS=0.000000;SumGene=0.000000;rpkmTSS=0.000000;rpkmGene=0.000000
Chr1    pausing_index   gene    712474  726891  0.788316        -       .       Parent=AT1G03060.1;SumExonLenth=11304;SumTSS=14.000000;SumGene=986.000000;rpkmTSS=5.263243;rpkmGene=6.676561
Chr1    pausing_index   gene    729935  731557  14.4008 +       .       Parent=AT1G03070.1;SumExonLenth=1019;SumTSS=211.000000;SumGene=60.000000;rpkmTSS=79.324591;rpkmGene=5.508365
Chr1    pausing_index   gene    729901  731525  12.5728 +       .       Parent=AT1G03070.2;SumExonLenth=1054;SumTSS=212.000000;SumGene=72.000000;rpkmTSS=79.700537;rpkmGene=6.339135
Chr1    pausing_index   gene    731612  738056  .       -       .       Parent=AT1G03080.1;SumExonLenth=5736;SumTSS=0.000000;SumGene=0.000000;rpkmTSS=0.000000;rpkmGene=0.000000
Chr1    pausing_index   gene    739679  744184  11.73   +       .       Parent=AT1G03090.1;SumExonLenth=2546;SumTSS=3.000000;SumGene=3.000000;rpkmTSS=1.127838;rpkmGene=0.096150
Chr1    pausing_index   gene    739679  744184  12.03   +       .       Parent=AT1G03090.2;SumExonLenth=2606;SumTSS=3.000000;SumGene=3.000000;rpkmTSS=1.127838;rpkmGene=0.093752
Chr1    pausing_index   gene    743885  746479  .       -       .       Parent=AT1G03100.1;SumExonLenth=2595;SumTSS=0.000000;SumGene=0.000000;rpkmTSS=0.000000;rpkmGene=0.000000
Chr1    pausing_index   gene    747197  748057  .       +       .       Parent=AT1G03103.1;SumExonLenth=640;SumTSS=0.000000;SumGene=0.000000;rpkmTSS=0.000000;rpkmGene=0.000000

```

In the output file the reference name will be compatible with the annotation. You may wish to change this later using `sed`.

### Additional info

Notice that once complete it will give you the sum of the score from the input w1 GFF. This is so that you can check any calculations by hand.

```
$ ./pausing_index -i example_data/example.w1.gff -a example_data/example_annotation.gff -tss 200 -o out.gff
The sum of the score (column 6) from example_data/example.w1.gff was 869858.21042313427
```

## Calculation

```
tss_rpkm = ((score_for_tss*1000000000)/sum_score) / tss_size;
gene_rpkm = ((score_for_gene*1000000000)/sum_score) / (sum_of_exons - tss_size);
```

## Install

Getting this to compile on work machine;

```
~/cmake-3.6.0-rc2-Linux-x86_64/bin/cmake ../pausing_index -DCMAKE_MODULE_PATH=../../seqan/util/cmake -DSEQAN_INCLUDE_PATH=../../seqan/include -DCMAKE_CXX_FLAGS=-std=c++14 -DCMAKE_BUILD_TYPE=Release
```

