#!/bin/bash
cd /project/functional-genomics/2019/group1/mapping/48hMEF/ChipSeq/
 for filename in *_filtered.sorted.bam
 do
   first=$(echo "$filename")
   first2=${first%_*}
   echo "$first"
   echo "$filename"
   cd /package/sequencer/anaconda3/envs/ngs/bin/
   ./plotFingerprint -b /project/functional-genomics/2019/group1/mapping/48hMEF/ChipSeq/"$filename" /project/functional-genomics/2019/group1/mapping/48hMEF/ChipSeq/48h_Inputnative_MNase_ChIP-Seq_filtered.sorted.bam \
   --labels "$first2" control \
   -plot /project/functional-genomics/2019/group1/fingerprints/48h/"$first2".png \
   --outRawCounts /project/functional-genomics/2019/group1/fingerprints/48h/"$first2".tab \
   --outQualityMetrics /project/functional-genomics/2019/group1/fingerprints/48h/"$first2".txt \
   --numberOfProcessors 3
done
