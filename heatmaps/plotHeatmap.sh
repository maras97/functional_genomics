#!/bin/bash

cd /project/functional-genomics/2019/group1/heatmap_data/48h

for filename in *matrix.gz
do
	mn=$(echo "$filename")
	mn=${mn:0:9}
	echo "$mn"
	/package/sequencer/anaconda3/envs/ngs/bin/plotHeatmap -m "$filename" -out "heatmap_""$mn"".png" --sortRegions descend
done


