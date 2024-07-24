#!/bin/bash
# This script generates and filters bam files. 
# Scroll to bottom of script for more information on SAM/BAM format, more particulary sam flag, alignment score and mapping quality

# Note, this script does not apply indel realignment.
# It is assumed that genotype calling from the final bam-files will be performed using bcftools mpileup (with the BAQ score calculating not disabled. I.e. do NOT set the option -B, --no-BAQ.)

#######################
# User-defined section

# PATHS TO FILES AND SOFTWARE:
REFERENCE=/home/mdejong/bearproject/refgenome/brownbear_chrom/ASM358476v1_HiC.fasta
REFBN=/home/mdejong/bearproject/refgenome/brownbear_chrom/ASM358476v1_HiC

REPEATFILE=myrepeats.bed                	# file needed for optional repeatmasking. This bed-file can be generated by editing output of repeat masking software.

INPUTDIR=/home/mdejong/bearproject/QCout
OUTPUTDIR=/home/mdejong/bearproject/bwaoutput_2024

BWA=/opt/software/bwa/bwa
SAMTOOLS=/opt/software/samtools-1.20/bin/samtools
PICARD=/opt/software/picard/picard.jar
BEDTOOLS=bedtools

NRCORES=14					# Number of cores or threads per file. For example: if working on 15 samples, and if specifying 4 cores, then usage of up to 60 cores.

# CONTROL PANEL: 				# Set the following flags one by one to TRUE. Not all the same time!
indexref=FALSE					# step 0. index reference
extractbn=FALSE					# step 0. extract base names of input files
mapdata=FALSE					# step 1. map with bwa and sort with samtools
getstats=FALSE					# step 2. get mapping statistics
indexbam=FALSE					# step 3. index bam files
addRG=FALSE					# step 4. add readgroups	
removedup=FALSE					# step 5: mark duplicates with PICARD markduplicates
filterbam=FALSE					# step 6: filter on mapping quality and alignment score
indexbam2=TRUE					# step 7: index filtered bam files
removerepeats=FALSE				# step 8: remove reads which mapped to a repetitive region
indexbam3=FALSE					# step 9: index repeatmasked bam files
#######################


# OVERVIEW OF STEPS:
# step 0. index reference genome 					
# wget link 
# $BWA index -a bwtsw $REFERENCE

# step 1a. map with bwa							
# step 1b. sort bam with samtools
# $BWA mem $REFERENCE ${INPUTDIR}/${bn1}_1.fq.gz ${INPUTDIR}/${bn1}_2.fq.gz -t 4 | samtools sort -o ${OUTPUTDIR}/${bn2}.sorted.bam	

# step 2. index bam with samtools (and get mapping statistics)		
# $SAMTOOLS index ${OUTPUTDIR}/${bn}.sorted.bam

# step 4. add readgroups											
# java -Xmx500g -jar $PICARD AddOrReplaceReadGroups I=${OUTPUTDIR}/${bn}.sorted.bam O=${OUTPUTDIR}/${bn}.sorted.RG.bam RGID=${bn} RGPL=illumina RGPU=unit1 RGLB=${bn} RGSM=${bn}

# step 5: mark duplicates with PICARD markduplicates			
# java -jar $PICARD MarkDuplicates I=${OUTPUTDIR}/${bn}.sorted.RG.bam O=${OUTPUTDIR}/${bn}.sorted.RG.dupremoved.bam M=${OUTPUTDIR}/mark_dup_metrics.txt REMOVE_DUPLICATES=True

# step 6. filter bam on mapping quality and/alignment score
# only needed if you will not set a mapping quality filter during genotype calling, or if you want to filter on additional settings:
# $SAMTOOLS view -b -h -q 20 -e '[AS]>=100' ${OUTPUTDIR}/${bn}.sorted.RG.dupremoved.bam > ${OUTPUTDIR}/${bn}.sorted.RG.dupremoved.filtered.bam

# step 7. index filtered bam files						
# $SAMTOOLS index ${OUTPUTDIR}/${bn}.sorted.RG.dupremoved.filtered.bam

# step 8. remove repetitive regions (optional)					
# $BEDTOOLS intersect -abam  ${INPUTDIR}/${bn}.sorted.RG.dupremoved.filtered.bam -b $REPEATFILE -v > ${OUTPUTDIR}/${bn}.sorted.RG.dupremoved.filtered.repeatmasked.bam &

# step 9. index filtered bam files						
# $SAMTOOLS index ${OUTPUTDIR}/${bn}.sorted.RG.dupremoved.filtered.repeatmasked.bam

# optional step. quality control						
# qualimap bamqc -bam $file -outdir ${OUTPUTDIR}/QC_${bn} -outfile ${bn}_qualimap2.pdf -nt 15 --java-mem-size=15G

# optional step. create consensus fasta
# angsd -i ${bn}.sorted.RG.markdup.realigned.bam -doFasta 3 -doCounts 1 -out ${bn}.fasta




if [[ "$indexref" = TRUE ]]
	then
	echo "Indexing reference..."
	bwa index -a bwtsw $REFERENCE				# if you map with bwa
	# bowtie2-build ${REFBN}.fna ${REFBN}		# if you map with bowtie2
	wait
	echo "Finished indexing."
fi

if [[ "$extractbn" = TRUE ]]
	then
	ls -1 ${INPUTDIR}/*_1.fq.gz > ${OUTPUTDIR}/allsamples.txt
	rev ${OUTPUTDIR}/allsamples.txt | cut -d "/" -f1 | rev | sed 's/.QC_1.fq.gz//g' > ${OUTPUTDIR}/allsamples.bn.txt
	rev ${OUTPUTDIR}/allsamples.txt | cut -d "/" -f1 | rev | sed 's/.QC_1.fq.gz//g' > allsamples.bn.txt
	echo "Created in input and output directory a file called 'allsamples.bn.txt', which lists the basenames of all input files."
	#
	# Optionally split in smaller files to avoid overloading system:
	split -l 8 ${OUTPUTDIR}/allsamples.bn.txt ${OUTPUTDIR}/subset
	split -l 8 ${OUTPUTDIR}/allsamples.bn.txt subset
	ls -1 ${OUTPUTDIR}/subset* > ${OUTPUTDIR}/mysubsets.txt
	echo "Created in input and output directory files called 'subset', with subsets of input files."
	fi

if [[ "$mapdata" = TRUE ]]
	then
	echo "Mapping with bwa..."
	for bn in $(cat allsamples.bn.txt)
		do
		echo $bn
		$BWA mem $REFERENCE ${INPUTDIR}/${bn}.QC_1.fq.gz ${INPUTDIR}/${bn}.QC_2.fq.gz -t $NRCORES | ${SAMTOOLS} sort --threads 3 -o ${OUTPUTDIR}/${bn}.sorted.bam &
		#$BWA mem $REFERENCE ${INPUTDIR}/${bn}_1.fq.gz ${INPUTDIR}/${bn}_2.fq.gz -t $NRCORES | ${SAMTOOLS} sort --threads 5 -o ${OUTPUTDIR}/${bn}.sorted.bam
		# -M: set shorter split reads as secondary. Might be better to include, but since I did not do it for most deer and bear samples, for consistency I will not include this option.  
		# or:
		# echo "Mapping with bowtie..."
		# bowtie2 -t -p 15 -q -x ${REFBN} -1 ${INPUTDIR}/${bn}.filteredreads_1P.fq.gz -2 ${INPUTDIR}/${bn}.filteredreads_2P.fq.gz -S ${OUTPUTDIR}/${bn}.sam &
		# ${SAMTOOLS} sort ${OUTPUTDIR}/${bn}.bam > ${OUTPUTDIR}/${bn}.sorted.bam
		done
	wait
	echo "Finished mapping data."
fi
	
if [[ "$getstats" = TRUE ]]
	then
	echo "Generating mapping statistics..."
	for bn in $(cat allsamples.bn.txt)
		do
		echo $bn
		${SAMTOOLS} flagstat ${OUTPUTDIR}/${bn}.sorted.bam > ${OUTPUTDIR}/${bn}.mappingstats.txt &
		done
	wait
	echo "Mapping statistics have been saved in files ending on 'mappingstats.txt'."
fi

if [[ "$indexbam" = TRUE ]]
	then
	echo "Indexing bam..."
	for bn in $(cat allsamples.bn.txt)
		do
		echo $bn
		${SAMTOOLS} index ${OUTPUTDIR}/${bn}.sorted.bam -@ 4 &
		# @: number of threads (although rarely above 130% anyway)
		done
	wait
	echo "Finished indexing."
fi	

if [[ "$addRG" = TRUE ]]
	then
	# uses up to 130% per sample, and 0.1% memory. So many samples can be run at the same time.
	# Around 2-3 hours per sample
	echo "Adding read groups..."
	for bn in $(cat allsamples.bn.txt)
		do
		echo $bn
		java -jar ${PICARD} AddOrReplaceReadGroups I=${OUTPUTDIR}/${bn}.sorted.bam O=${OUTPUTDIR}/${bn}.sorted.RG.bam RGID=${bn} RGPL=illumina RGPU=unit1 RGLB=${bn} RGSM=${bn} &
		# we set RGPL to illumina because novogene uses the illumina platform. we don't define platform unit (RGPU) 
		# rm ${OUTPUTDIR}/${bn}.sorted.bam 
		# rm ${OUTPUTDIR}/${bn}.sorted.bam.bai
		done
	wait
	echo "Finished adding read groups..."
fi

if [[ "$removedup" = TRUE ]]
	then
	echo "Marking and removing duplicates..."
	# next command spikes now and then up to 50 cores per sample (if available). 
	# Plus, depending on input bam-file size, can consume a lot of memory (in Cetus around 15% per sample)  
	# Therefore, I run maximum 3 samples at a time. 
	# the argument MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 is included to avoid the error: 
	# Exception in thread "main" htsjdk.samtools.SAMException: /tmp/mdejong/CSPI.6051568413464160692.tmp/18682.tmpnot found. Caused by: java.io.FileNotFoundException: /tmp/mdejong/CSPI.6051568413464160692.tmp/18682.tmp (Too many open files)
	for bn in $(cat allsamples.bn.txt)
		do
		echo $bn
		java -Xmx64g -Djava.io.tmpdir=picardtmp -jar ${PICARD} MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 I=${OUTPUTDIR}/${bn}.sorted.RG.bam O=${OUTPUTDIR}/${bn}.sorted.RG.dupremoved.bam M=${OUTPUTDIR}/${bn}.mark_dup_metrics.txt REMOVE_DUPLICATES=True
		#rm ${OUTPUTDIR}/${bn}.sorted.RG.bam
		done
	wait
	echo "Finished marking and removing duplicates."
fi

if [[ "$filterbam" = TRUE ]]
	then
	# takes approximately 1 hour per sample; assuming 3 cores can do around 20 samples at a time without clogging the system
	# q:		filter on mapping quality (--min-MQ)
	# h:		include header
	# b:		output in bam format
	# e:		filter on expression; '[AS]>=100' means: keep reads with an alignment score above 100		
	# f 0x2:	keep properly paired reads only
	# F4:		remove unmapped read
	# [AS]>=100	keep only reads with an alignment score of 100 or higher.  
	#		Alignment score is calculated as follows: a match is 1, a mismatch is -4, a gap opening is -6, and a gap extension is -1. 
	#		Also note that alignment score considers base qualities (if I am not mistaken).
	#		So you also can get a penalty if base quality score is below 13 (default value), meaning:  !"#$%&'()*+,-. (not sure about the period)
	echo "Filtering on mapping quality and alignment score..."
	#for bn in $(cat subsetab)
	for bn in $(cat allsamples.bn.txt)
		do
		echo $bn
		${SAMTOOLS} view -@ $NRCORES -bhq 20 -F4 -f 0x2 -e '[AS]>=100' ${OUTPUTDIR}/${bn}.sorted.RG.dupremoved.bam > ${OUTPUTDIR}/${bn}.sorted.RG.dupremoved.filtered.bam &
                #${SAMTOOLS} view -@ $NRCORES -bhq 20 -F4 -e '[AS]>=35' ${OUTPUTDIR}/${bn}.sorted.RG.dupremoved.bam > ${OUTPUTDIR}/${bn}.sorted.RG.dupremoved.filtered.AS35.bam &
		#rm ${OUTPUTDIR}/${bn}.sorted.RG.dupremoved.bam
		done
	wait
	echo "Finished filtering bam files."
fi

if [[ "$indexbam2" = TRUE ]]
	then
	echo "Indexing bam file..."
	for bn in $(cat allsamples.bn.txt)
		do
		echo $bn
		${SAMTOOLS} index ${OUTPUTDIR}/${bn}.sorted.RG.dupremoved.filtered.bam -@ 3 &
		done
	wait
	echo "Finished indexing bam file."
fi

if [[ "$removerepeats" = TRUE ]]
        then
        echo "Removing reads which mapped to repetitive regions..."
        for bn in $(cat allsamples.bn.txt)
            do
            echo $bn
            ${BEDTOOLS} intersect -abam  ${INPUTDIR}/${bn}.sorted.RG.dupremoved.filtered.bam -b $REPEATFILE -v > ${OUTPUTDIR}/${bn}.sorted.RG.dupremoved.filtered.repeatmasked.bam &
            done
        wait
        echo "Reads have been removed. New output files have the suffix: 'sorted.RG.dupremoved.filtered.repeatmasked.bam'."
fi

if [[ "$indexbam3" = TRUE ]]
	then
	echo "Indexing bam file..."
	for bn in $(cat allsamples.bn.txt)
		do
		echo $bn
		${SAMTOOLS} index ${OUTPUTDIR}/${bn}.sorted.RG.dupremoved.filtered.repeatmasked.bam -@ 3 &
		done
	wait
	echo "Finished indexing bam file."
fi


# OPTIONALLY GENERATING QC REPORTS:

# echo "Validating bam file"
# java -jar ${PICARD} ValidateSamFile I=${OUTPUTDIR}/${bn}.sorted.RG.dupremoved.filtered.bam MODE=SUMMARY

# echo "Generating QC report..."
# Optionally install qualimap:
# conda activate conenv3.7
# conda install -c bioconda qualimap
	
# mkdir ${OUTPUTDIR}/QC_${bn}
# qualimap bamqc -bam $file -outdir ${OUTPUTDIR}/QC_${bn} -outfile ${OUTPUTDIR}/${bn}_qualimap2.pdf -nt 10 --java-mem-size=10G








############## UNDERSTANDING SAM/BAM-FORMAT ###############

# SAM FLAG (2nd column)
# Binary number composed of following elements:
# 0x1:		1		read paired
# 0x2: 		2		read mapped in proper pair
# 0x4:		4		read unmapped
# 0x8:		8		mate unmapped
# 0x16:		16		read reverse strand
# 0x20:		32		mate reverse strand
# 0x40: 	64		first in pair
# 0x80: 	128		second in pair
# 0x100: 	256		not primary alignment
# 0x200: 	512		read fails platform/vendor quality checks
# 0x400: 	1024	read is PCR or optical duplicates
# 0x800:	2048	supplementary alignment
 

# MAPPING QUALITY (5th column)
# Phred-scaled probability that read has been mapped to the wrong position in genome:
# MAPQ = -10*log10(pr_wrongly_mapped)
# Roughly: 
# Q10 = Pr(correctly mapped) = 90%. 
# Q20 = Pr(correctly mapped) = 99%. 
# Q30 = Pr(correctly mapped) = 99.9%. 
# Etc.
# BWA: 0-60; Bowtie: 0-42.
# In theory, probability depends on:
# 1. Contamination/sequencing errors
# 2. Mapping algorithm heuristics
# 3. Error due to the repetiveness of reference
# The MAPQ only considers the third factor.
# It compares the alignment score of the best alignment (AS:i) to the alignment score of the second-best alignment (XS:i). (See below for more info on AS and XS).
# If AS is much higher than XS, low probability of wrong mapping position (and thus high MAPQ value, close to 60).
# If AS==XS, MAPQ = 0.
# Note the MAPQ-score does not tell anything about the alignment itself.
# You can have a very shitty alignment (many gaps, many mismatches), but if the second best alignment is even much more shitty, then still high MAPQ score.
# If there is no second mapping position, MAPQ is set to maximum (60) in the case of BWA, and 255 in the case of Bowtie (?). 

# CIGAR STRING (6th column)
# describes alignment. Two examples:
# 150M:				optimal score for a read of 150bp, namely 150 matches/mismatches
# 23M2D50M2I57M6S:	23 matches/mismatches, 2 deletions, 50 matches, 2 insertions, 59 matches, 18 soft-clipped sites (segment at end of read that does not appear in alignment)  	

# AS:i:count
# alignment score.
# In case of BWA very simple: simply the number of matches. In the examples for CIGAR above, the alignment scores are respectively 150 and 132 (23 + 50 + 59).
# Note that 'i' simply indicates that this is an integer, and not a code (Z) for example.
 
# XS:i:count
# alignment score of second best alignment, if present

# MC:Z
# CIGAR string of second best alignment (or read pair mate?)

# MD:Z
# CIGAR string for reference (rather than query)

# NM:i:count
# Number of mismatches/deletions/insertions

# NH:i:2
# Number of regions in reference to which query could be aligned to








