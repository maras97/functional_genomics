#!/bin/bash
cd /project/functional-genomics/2019/data/sra/48hMEF/RNASeq/
for file in *.fastq
do
  first=$(echo "$file")
  first=${first:0:10}
  cd /project/functional-genomics/2019/group1/mapping/48hMEF/RNASeq/
  STAR --runThreadN 10 \
  --runMode alignReads \
  --genomeDir /project/functional-genomics/2019/data/genome/STARindex/ \
  --readFilesIn /project/functional-genomics/2019/data/sra/48hMEF/RNASeq/"$file" \
  --sjdbGTFfile  /project/functional-genomics/2019/data/annotation/Mus_musculus.NCBIM37.67.gtf \
  --outSAMtype BAM SortedByCoordinate \
  --outWigType wiggle \
  --outWigNorm RPM \
  --outFileNamePrefix ./"$first" \
  --quantMode GeneCounts
done
