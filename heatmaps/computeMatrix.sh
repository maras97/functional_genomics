#!/bin/bash

cd /project/functional-genomics/2019/group1/heatmap_data/48h

for filename in *intersect.bed
do
	fn=$(echo "$filename")
	f1=${fn:0:4}
	f2=${fn:5:4}
	/package/sequencer/anaconda3/envs/ngs/bin/computeMatrix reference-point --referencePoint center -b 1000 -a 1000 -R "$f1""_wo_""$f2"".bed" "$f2""_wo_""$f1"".bed" "$filename" -S "48h_""$f1"".bw" "48h_""$f2"".bw" -p 15 -o "$f1""_""$f2""_matrix.gz"
done
