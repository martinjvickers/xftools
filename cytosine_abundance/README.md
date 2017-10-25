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

## Downstream analysis

The output looks like this;

```
1       xftools window  50      100     16      .       .       CAA=0;CAC=0;CAG=0;CAT=3;CCA=1;CCC=3;CCG=0;CCT=4;CG=0;CGA=0;CGC=0;CGG=0;CGN=0;CGT=0;CHG=0;CHH=16;CNG=0;CNN=0;CTA=4;CTC=0;CTG=0;CTT=1
```

So if you want to see just the score and the three contexts (CG, CHG and CHH) you can do something like;

cat meh.gff  | awk '{split($9,a,";"); print $6"\t"a[9]"\t"a[15]"\t"a[16]}'

### Context layout

|Column | Context  |
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


`
CAA=0
CAC=0
CAG=0
CAT=0
CCA=0
CCC=0
CCG=0
CCT=0
CG=0
CGA=0
CGC=0
CGG=0
CGN=0
CGT=0
CHG=0
CHH=0
CNG=0
CNN=0
CTA=0
CTC=0
CTG=0
CTT=0
```
