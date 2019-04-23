#!/bin/bash
cd /project/functional-genomics/2019/data/sra/48hMEF/ATACSeq
for filename in *.fastq;
do

  echo "$filename"
  first=$(echo "$filename" )
  first=${first:0:10}
  echo "$first"

  cd /project/functional-genomics/2019/data/genome

  touch /project/functional-genomics/2019/group1/mapping/48hMEF/ATACSeq/"$first".sam


  bowtie2 -x mm9 -t -p 20 -U /project/functional-genomics/2019/data/sra/48hMEF/ATACSeq/"$filename" -S /project/functional-genomics/2019/group1/mapping/48hMEF/ATACSeq/"$first".sam

  cd /project/functional-genomics/2019/group1/mapping/48hMEF/ATACSeq/

  samtools view -@16 -bS "$first".sam > "$first".bam

  samtools sort -@16 "$first".bam -o "$first".sorted.bam

  samtools flagstat -@16 "$first".sorted.bam > "$first"_flagstat.txt


done
