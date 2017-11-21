# methylation_tools

## trim

### Compiling

Getting this to compile on work machine;

```
~/cmake-3.6.0-rc2-Linux-x86_64/bin/cmake ../trim -DCMAKE_MODULE_PATH=../seqan/util/cmake -DSEQAN_INCLUDE_PATH=../seqan/include -DCMAKE_CXX_FLAGS=-std=c++14 -DCMAKE_BUILD_TYPE=Release
```

```
mvickers@n108379:~/development/xftools/trim$ bowtie2 -k 4 -X 4000 --local --threads 2 -x TAIR10_reference_wCM.fa -1 H1VN_rep1_R1_trim.fastq.gz -2 H1VN_rep1_R2_trim.fastq.gz -S H1VN_rep1.sam
Warning: gzbuffer added in zlib v1.2.3.5. Unable to change buffer size from default of 8192.
Warning: gzbuffer added in zlib v1.2.3.5. Unable to change buffer size from default of 8192.
2557801 reads; of these:
  2557801 (100.00%) were paired; of these:
    258665 (10.11%) aligned concordantly 0 times
    1548094 (60.52%) aligned concordantly exactly 1 time
    751042 (29.36%) aligned concordantly >1 times
    ----
    258665 pairs aligned concordantly 0 times; of these:
      81054 (31.34%) aligned discordantly 1 time
    ----
    177611 pairs aligned 0 times concordantly or discordantly; of these:
      355222 mates make up the pairs; of these:
        287161 (80.84%) aligned 0 times
        10569 (2.98%) aligned exactly 1 time
        57492 (16.18%) aligned >1 times
94.39% overall alignment rate
mvickers@n108379:~/development/xftools/trim$ samtools view -bS H1VN_rep1.sam > H1VN_rep1.bam
mvickers@n108379:~/development/xftools/trim$ samtools sort H1VN_rep1.bam > H1VN_rep1_sorted.bam
mvickers@n108379:~/development/xftools/trim$ samtools index H1VN_rep1_sorted.bam

```

```
# =============================
# Remove  unmapped, mate unmapped
# not primary alignment, reads failing platform
# ==================  
mvickers@n108379:~/development/xftools/trim$ samtools view -h H1VN_rep1_sorted.bam | ./assign_multimappers.py -k 4 | samtools view -h -F 1804 | samtools view -bS > H1VN_rep1_filtered.bam
```

```
# ========================
# Mark duplicates
# ======================
[mvickers@NBI-HPC interactive for_shengbo]$ picard MarkDuplicates INPUT=H1VN_rep1_filtered.bam OUTPUT=H1VN_rep1_filtered_dupmark.bam METRICS_FILE=H1VN_rep1_filtered_dup.qc VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false
[Mon Jul 10 14:29:04 BST 2017] picard.sam.markduplicates.MarkDuplicates INPUT=[H1VN_rep1_filtered.bam] OUTPUT=H1VN_rep1_filtered_dupmark.bam METRICS_FILE=H1VN_rep1_filtered_dup.qc REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT    MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 SORTING_COLLECTION_SIZE_RATIO=0.25 PROGRAM_RECORD_ID=MarkDuplicates PROGRAM_GROUP_NAME=MarkDuplicates DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES READ_NAME_REGEX=[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).* OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 VERBOSITY=INFO QUIET=false COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json
[Mon Jul 10 14:29:04 BST 2017] Executing as mvickers@n128n1.nbicluster on Linux 3.10.0-514.26.1.el7.x86_64 amd64; Java HotSpot(TM) 64-Bit Server VM 1.8.0_45-b14; Picard version: 1.134(a7a08c474e4d99346eec7a9956a8fe71943b5d80_1434033355) IntelDeflater
INFO	2017-07-10 14:29:05	MarkDuplicates	Start of doWork freeMemory: 1999806904; totalMemory: 2022178816; maxMemory: 7635730432
INFO	2017-07-10 14:29:05	MarkDuplicates	Reading input file and constructing read end information.
INFO	2017-07-10 14:29:05	MarkDuplicates	Will retain up to 29368193 data points before spilling to disk.
INFO	2017-07-10 14:29:25	MarkDuplicates	Read     1,000,000 records.  Elapsed time: 00:00:18s.  Time for last 1,000,000:   18s.  Last read position: Chr1:28,444,702
INFO	2017-07-10 14:29:25	MarkDuplicates	Tracking 23772 as yet unmatched pairs. 21969 records in RAM.
INFO	2017-07-10 14:29:36	MarkDuplicates	Read     2,000,000 records.  Elapsed time: 00:00:29s.  Time for last 1,000,000:   10s.  Last read position: Chr3:4,697,904
INFO	2017-07-10 14:29:36	MarkDuplicates	Tracking 90462 as yet unmatched pairs. 3588 records in RAM.
/nbi/software/testing/bin/core/../..//picard/1.134/x86_64/bin/picard: line 3: 42833 Killed                  java -Xmx8g -jar /nbi/software/testing/picard/1.134/x86_64/jars/picard.jar $@
```

