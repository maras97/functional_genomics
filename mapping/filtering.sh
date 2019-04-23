#!/bin/bash
cd /project/functional-genomics/2019/group1/mapping/48hMEF/ChipSeq

for file in *_filtered.sorted.bam
do
    echo "$file"
    first=$(echo "$file")
    first=${first%_*}
    samtools rmdup -S /project/functional-genomics/2019/group1/mapping/48hMEF/ChipSeq/"$file" /project/functional-genomics/2019/group1/mapping/48hMEF/Filtered_ChipSeq/"$first".bam
done
