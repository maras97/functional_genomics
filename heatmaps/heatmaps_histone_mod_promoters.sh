#!/bin/bash

/package/sequencer/anaconda3/envs/ngs/bin/computeMatrix reference-point --referencePoint center -b 2000 -a 2000 \
      -R /project/functional-genomics/2019/group1/up_ESC.bed /project/functional-genomics/2019/group1/down_ESC.bed  -S /project/functional-genomics/2019/group1/mapping/ESC/Filtered_ChipSeq/temp/*.bw \
      -p 20 -o /project/functional-genomics/2019/group1/HM_HighLowESC.gz

/package/sequencer/anaconda3/envs/ngs/bin/plotHeatmap -m /project/functional-genomics/2019/group1/HM_HighLowESC.gz \
        --samplesLabel "H3.3" "H3K27ac" "H3K27me3" "H3K36me3" "H3K4me1" "H3K4me2" "H3K4me3" "H3K79me2" "H3K9ac" "H3K9me3" "H3" "Hdac1"\
        --regionsLabel "up" "down"\
        -out /project/functional-genomics/2019/group1/fotos/Heatmaps/heatmap_HM_HighLowESC.png \
        --sortRegions descend
