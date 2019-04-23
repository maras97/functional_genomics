#!/bin/bash

cd /project/functional-genomics/2019/group1/mapping/ESC/ChipSeq/

for filename in ESC_H3K4me2_ChIP-Seq_filtered.sorted.bam
do
    first=$(echo "$filename")
    first=${first%_filtered*}

    cd /project/functional-genomics/2019/data/genome
    macs2 callpeak -t /project/functional-genomics/2019/group1/mapping/ESC/ChipSeq/"$filename" -c /project/functional-genomics/2019/group1/mapping/ESC/ChipSeq/ESC_Inputnative_MNase_ChIP-Seq_filtered.sorted.bam -f BAM -g mm --bw 150 -q 0.005 --mfold 4 50 --keep-dup 1 --outdir /project/functional-genomics/2019/group1/peak_calling/ESC/ChipSeq/ -n "$first"
done
