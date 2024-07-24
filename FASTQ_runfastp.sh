#!/bin/bash

########
FASTP=/home/mdejong/software/fastp
QCDIR=/home/mdejong/bearproject/QCout
########

for file in $(cat allsamples.bn.txt)
do
bn=$(echo $file | sed 's/_1.fq.gz//g')
echo $bn
$FASTP --in1 ${bn}_1.fq.gz --in2 ${bn}_2.fq.gz --out1 ${QCDIR}/${bn}.QC_1.fq.gz --out2 ${QCDIR}/${bn}.QC_2.fq.gz --thread 14 --detect_adapter_for_pe --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 --qualified_quality_phred 15 --unqualified_percent_limit 40 --n_base_limit 5 --length_required 36 --low_complexity_filter --correction --overrepresentation_analysis 2> ${QCDIR}/${bn}.fastp.log.txt &
done
wait
echo "Analysis finished."
