# XFTOOLS

## TOOLNAME

Description of at_compare. Given two GFF files and a reference genome,
this program will output the GFFs with respective feature AT value as a tag
appended to the end of the GFF. e.g.

```
chr1	TAIR10	TE_Tair10	11897	11976	.	+	.	ID=AT1TE00010;Name=AT1TE00010;Alias=ATCOPIA24;AT=0.6375
chr1	TAIR10	TE_Tair10	16883	17009	.	-	.	ID=AT1TE00020;Name=AT1TE00020;Alias=ATREP4;AT=0.669291
chr1	TAIR10	TE_Tair10	17024	18924	.	+	.	ID=AT1TE00025;Name=AT1TE00025;Alias=ATREP3;AT=0.776433
chr1	TAIR10	TE_Tair10	18331	18642	.	-	.	ID=AT1TE00030;Name=AT1TE00030;Alias=ATHATN7;AT=0.682692
```

### Compiling

Getting this to compile on work machine;

```
~/cmake-3.6.0-rc2-Linux-x86_64/bin/cmake ../at_compare \
   -DCMAKE_MODULE_PATH=../../seqan/util/cmake \
   -DSEQAN_INCLUDE_PATH=../../seqan/include \
   -DCMAKE_CXX_FLAGS=-std=c++14 \
   -DCMAKE_BUILD_TYPE=Release \
   -DCMAKE_CXX_ARGS="-DSEQAN_DISABLE_VERSION_CHECK=YES" 
```

### Running

```
./at_compare -i example_data/TAIR10_reference.fas \
-1 example_data/TAIR10_GFF3-gene_only_MV.gff \
-2 example_data/TAIR10_TE_Anno.XF.gff \
-o1 example_data/TAIR10_gene_at.gff \
-o2 example_data/TAIR10_TE_at.gff
```
