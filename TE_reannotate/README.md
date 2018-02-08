# XFTOOLS

## TE_reannotate

Description of TE_reannotate

So, taking the single-c GFF file, first we make continous blocks in the 
genome where methylation is > a variable. The default will be >0.0. 
Don't take into account strandedness at this point. 

Then, find all of the reads that overlap the block. This also includes 
the scenario where a read overlaps that block, and a secondary read 
doesn't overlap the block but overlaps the previous read. e.g.

       [---BLOCK---]  
   <---->
<---->

Finally, within the selected reads for that block, we group together reads 
into units that are within 200bp (make adjustable).

e.g.

       [-----------------BLOCK-------------]
     <----->                     <---->
   <---->                                 <---->
<---->                                        <---->

<--UNIT---->                     <------UNIT------->



### Compiling

Getting this to compile on work machine;

```
~/cmake-3.6.0-rc2-Linux-x86_64/bin/cmake ../TE_reannotate \
   -DCMAKE_MODULE_PATH=../../seqan/util/cmake \
   -DSEQAN_INCLUDE_PATH=../../seqan/include \
   -DCMAKE_CXX_FLAGS=-std=c++14 \
   -DCMAKE_BUILD_TYPE=Release \
   -DCMAKE_CXX_ARGS="-DSEQAN_DISABLE_VERSION_CHECK=YES" 
```



### Running


