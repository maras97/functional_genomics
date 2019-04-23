#!/bin/bash
cd /project/functional-genomics/2019/group1/mapping/48hMEF/ATACSeq
for filename in *.sorted.bam;
do
    first=$(echo "$filename" )
    first=${first:0:10}
    filtered=$(echo "$first""_filtered.sorted" )
    samtools view -bq 10 -@ 20 "$filename" > "$filtered".bam
    samtools index -@ 20 "$filtered".bam
    samtools idxstats -@ 20 "$filtered".bam > "$first".idxstats.txt
    samtools flagstat -@ 20 "$filtered".bam > "$first".flagstat.txt
    bamCoverage -b "$filtered".bam -o "$filtered".bw -of bigwig -bs 25 --normalizeUsing RPKM --ignoreForNormalization chrX--extendReads 200 --numberOfProcessors 3
done
                            
