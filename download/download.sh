#!/bin/bash

while IFS='' read -r line || [[ -n "$line" ]]; do
    prefetch $line
    fastq-dump --split-files --skip-technical  --readids --dumpbase --outdir /project/functional-genomics/2019/data/sra/ESC/ /project/functional-genomics/2019/data/sra/ESC/ncbi/public/sra/$line.sra

done < "$1"
