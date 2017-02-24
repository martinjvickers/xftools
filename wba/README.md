# methylation_tools

## Window By Annotation

### Compiling

Getting this to compile on work machine;

```
~/cmake-3.6.0-rc2-Linux-x86_64/bin/cmake ../wba -DCMAKE_MODULE_PATH=../../seqan/util/cmake -DSEQAN_INCLUDE_PATH=../../seqan/include -DCMAKE_CXX_FLAGS=-std=c++14 -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_ARGS="-DSEQAN_DISABLE_VERSION_CHECK=YES"
```

### Getting started


### some performance testing

####My run ~21mins

```
	Command being timed: "./wba -a example_data/captest.gff -i /mnt/group_share/Martin_Vickers/data_upload_for_papers/2016_SLM_paper/processed_data_to_upload/rdr2.spm.CHH.tair10.gff -o meh -l"
	User time (seconds): 1277.23
	System time (seconds): 0.77
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 21:19.06
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 5876
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 1497
	Voluntary context switches: 18
	Involuntary context switches: 204793
	Swaps: 0
	File system inputs: 0
	File system outputs: 2816
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0

```

####Original program ~35mins and died. I take this as a win for my program.

```
Command exited with non-zero status 255
	Command being timed: "/mnt/ssd/mvickers/rice_data/bin/window_by_annotation_mc.pl -g example_data/captest.gff -k -t ID -s -l -o meh.out /mnt/group_share/Martin_Vickers/data_upload_for_papers/2016_SLM_paper/processed_data_to_upload/rdr2.spm.CHH.tair10.gff"
	User time (seconds): 2155.13
	System time (seconds): 1.03
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 35:57.16
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 134048
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 33847
	Voluntary context switches: 21
	Involuntary context switches: 318838
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 255

```
