#!/bin/bash

cd /project/functional-genomics/2019/group1/heatmap_data/48h

for FN1 in *extd.bed;
do 
	for FN2 in *extd.bed;
	do
		bedtools intersect -a "$FN1" -b "$FN2" >> "${FN1:4:4}""_""${FN2:4:4}""_intersect.bed"
		bedtools intersect -a "$FN1" -b "$FN2" -v >> "${FN1:4:4}""_wo_""${FN2:4:4}"".bed"
		bedtools intersect -a "$FN2" -b "$FN1" -v >> "${FN2:4:4}""_wo_""${FN1:4:4}"".bed"
	done
done
