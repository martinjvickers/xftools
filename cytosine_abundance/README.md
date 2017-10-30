## cytosine_abundance

The purpose of this program is to calculate the cytosine abundance and cytosine context's in a reference genome.

## Usage

```
cytosine_abundance -i reference.fa -o reference_cytosine.gff -w 50
```

If you wish to use a GFF rather than a window size to calculate cytosine abundance then specify using the `-r` flag.

```
cytosine_abundance -i reference.fa -o reference_cytosine.gff -r region.gff
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
~/cmake-3.6.0-rc2-Linux-x86_64/bin/cmake ../cytosine_abundance -DCMAKE_MODULE_PATH=../seqan/util/cmake -DSEQAN_INCLUDE_PATH=../seqan/include -DCMAKE_CXX_FLAGS=-std=c++14 -DCMAKE_BUILD_TYPE=Release
```

## Downstream analysis

The output looks like this;

```
1       xftools window  50      100     16      .       .       CAA=0;CAC=0;CAG=0;CAT=3;CCA=1;CCC=3;CCG=0;CCT=4;CG=0;CGA=0;CGC=0;CGG=0;CGN=0;CGT=0;CHG=0;CHH=16;CNG=0;CNN=0;CTA=4;CTC=0;CTG=0;CTT=1
```

So if you want to see just the score and the three contexts (CG, CHG and CHH) you can do something like;

```
cat meh.gff  | awk '{split($9,a,";"); print $6"\t"a[9]"\t"a[15]"\t"a[16]}'
```

If you wish to see the score and the distribution of contexts for CG;

```
cat meh.gff  | awk '{split($9,a,";"); print $6"\t"a[9]"\t"a[10]"\t"a[11]"\t"a[12]"\t"a[13]"\t"a[14]}'
```

If you wish to see the score and the distribution of contexts for CHG;

```
cat meh.gff  | awk '{split($9,a,";"); print $6"\t"a[15]"\t"a[3]"\t"a[7]"\t"a[17]"\t"a[21]}'
```

If you wish to see the score and the distribution of contexts for CHH;

```
cat meh.gff  | awk '{split($9,a,";"); print $6"\t"a[16]"\t"a[1]"\t"a[2]"\t"a[4]"\t"a[5]"\t"a[6]"\t"a[8]"\t"a[18]"\t"a[19]"\t"a[20]"\t"a[22]"\t"}'
```

### Context position layout within column 9 of the GFF

Column 9 of the GFF output will adhere to this order.

|Position | Context  |
|-------|----------|
|1	|	CAA|
|2	|	CAC|
|3	|	CAG|
|4	|	CAT|
|5	|	CCA|
|6	|	CCC|
|7	|	CCG|
|8	|	CCT|
|9	|	CG|
|10	|	CGA|
|11	|	CGC|
|12	|	CGG|
|13	|	CGN|
|14	|	CGT|
|15	|	CHG|
|16	|	CHH|
|17	|	CNG|
|18	|	CNN|
|19	|	CTA|
|20	|	CTC|
|21	|	CTG|
|22	|	CTT|

