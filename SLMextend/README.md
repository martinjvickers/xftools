# XFTOOLS

## SLMextend

Description of SLMextend

### Compiling

Getting this to compile on work machine;

```
~/cmake-3.6.0-rc2-Linux-x86_64/bin/cmake ../SLMextend \
   -DCMAKE_MODULE_PATH=../seqan/util/cmake \
   -DSEQAN_INCLUDE_PATH=../seqan/include \
   -DCMAKE_CXX_FLAGS=-std=c++14 \
   -DCMAKE_BUILD_TYPE=Release \
   -DCMAKE_CXX_ARGS="-DSEQAN_DISABLE_VERSION_CHECK=YES" 
```

### Running

```
./SLMextend -slm example_data/converted/step5a_avg5_9612_SLM.gff -spm example_data/converted/sperm.CHH.w1.gff.gz -s1 example_data/converted/soma.9DAP_em.CHH.w1.gff.gz -s2 example_data/converted/SRR2079442_earshoot_rmdup_CHH.w1.gff.gz -s3 example_data/converted/SRR2079447_B73_shoot_apex_rmdup_CHH.w1.gff.gz
```
