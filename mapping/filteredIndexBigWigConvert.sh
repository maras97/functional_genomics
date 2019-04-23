#!/bin/bash
cd /project/functional-genomics/2019/group1/mapping/ESC/Filtered_ChipSeq/
for file in *.bam
do
    echo "$file"
    first=$(echo "$file")
    first=${first%.*}
    samtools index -@ 20 "$first".bam
    bamCoverage -b "$file" -o "$first".bw -of bigwig -bs 25 --normalizeUsing RPKM --ignoreForNormalization chrX --extendReads 200 --numberOfProcessors 20
done
