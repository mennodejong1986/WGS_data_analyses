#!/bin/bash
# This script is to extract regions from a bamfile
# the region(s) should be defined in a bed file

XBED=xchrom.bed
YBED=ychrom.bed

# ls -1 *sorted.RG.dupremoved.filtered.bam > mybamfiles.txt

for file in $(cat allsamples.bn.txt)
do
bn=$(echo $file | cut -f1 -d '.')
echo $bn
#samtools view -bh -L $XBED ${bn}.sorted.RG.dupremoved.filtered.bam -o ${bn}.xchrom.bam --threads 4 &
#samtools view -bh -L $YBED ${bn}.sorted.RG.dupremoved.filtered.bam -o ${bn}.ychrom.bam --threads 3 &
samtools index ${bn}.xchrom.bam &
samtools index ${bn}.ychrom.bam &
done
wait
echo "Finished subtracting region(s)."
