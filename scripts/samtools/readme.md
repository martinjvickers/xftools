
```
samtools view -b -f 16 accepted_hits.bam > reverse_strand.bam
samtools view -b -F 20 accepted_hits.bam > forward_strand.bam
```

##sort and view

```
samtools sort reverse_strand.bam > reverse_strand_sorted.bam
samtools index reverse_strand_sorted.bam
samtools sort forward_strand.bam > forward_strand_sorted.bam
samtools index forward_strand_sorted.bam
```

![alt text](images/igv_strands.png "Strands viewed in different tracks within IGV")

