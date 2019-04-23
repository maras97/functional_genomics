#!/bin/bash

cd /project/functional-genomics/2019/group1/peak_calling/48h/ChipSeq/

for filename in *.xls
do
	first=$(echo "$filename")
	clipped=${first%_peaks*}
	count=$(awk '($9>2.301){ ++count } END{print count}' $filename)
	avg=$(awk '{x+=$4; next} END{print x/NR}' $filename)
	echo -e "$clipped \t $count \t $avg" >> peak_count_width.txt	
done
