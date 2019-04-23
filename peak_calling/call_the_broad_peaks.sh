#!/bin/bash

cd /project/functional-genomics/2019/group1/mapping/ESC/ChipSeq/

for filename in ESC_H3K79me2_ChIP-Seq_filtered.sorted.bam ESC_H3K9me3_ChIP-Seq_filtered.sorted.bam ESC_H3K27me3_ChIP-Seq_filtered.sorted.bam ESC_H3K36me3_ChIP-Seq_filtered.sorted.bam
do
    first=$(echo "$filename")
    first=${first%_filtered*}	
    cd /project/functional-genomics/2019/data/genome
    macs2 callpeak -t /project/functional-genomics/2019/group1/mapping/ESC/ChipSeq/"$filename" -c /project/functional-genomics/2019/group1/mapping/ESC/ChipSeq/ESC_Inputnative_MNase_ChIP-Seq_filtered.sorted.bam -f BAM -g mm --bw 150 -q 0.005 --mfold 4 50 --keep-dup 1 --broad --outdir /project/functional-genomics/2019/group1/peak_calling/ESC/ChipSeq/broad_peaks/ -n "$first"
done
