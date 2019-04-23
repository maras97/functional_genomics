#!/bin/bash
cd /project/functional-genomics/2019/group1/mapping/48hMEF/RNASeq/
for file in *.wig
do
    name=$(echo "$file")
    name=${name%.wig}
    echo "$name"
    wigToBigWig "$file" /project/functional-genomics/2019/data/genome/STARindex/chrNameLength.txt "$name".bw
done
