#!/bin/bash

cd /project/functional-genomics/2019/group1/heatmap_data/48h

for FN in *.bed
do 	
	FN_NEW=$(echo "$FN")
	FN_NEW=${FN_NEW:0:8}
	FN_NEW=$(echo "$FN_NEW""_extd")
	awk '{print $1 "\t" ($2 - 100) "\t" ($3 + 100) "\t" $4 "\t" $5}' "$FN"  >> "$FN_NEW".bed	
done
