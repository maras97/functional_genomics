for file in /project/functional-genomics/2019/data/sra/48hMEF/*.fastq
do
  fastqc "$file"
done
