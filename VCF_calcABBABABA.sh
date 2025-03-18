#!/bin/bash
# script for sliding window ABBA-BABA analyses (D and f statistics)
# Using software developed by Simon Martin [Martin et al. 2015, MBE]

# Online tutorial:
# http://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/abba-baba-statistics/





### INSTALLATION
# Do not copy the executables of Simon Martin's software to your working directory, because they may depend on other executables. 
# For example, parseVCF.py depends on executable genomics.py. So keep all together in the git clone directory, using this command:
git clone https://github.com/simonhmartin/genomics_general




#### LOAD CONDA ENVIRONMENT
# Note: the tutorial mentioned that this software runs in a python 2.7 environment, but this seems to be incorrect and leads to errors.
# The software requires python3.7. 

conda create --name conenv3.7		# create new environment named ‘conenv3.7’
conda activate conenv3.7		# load environment
conda install numpy			# install software named numpy 




#### CONVERT DATA:
# Input data should be a gzipped vcf file with snp data.
python /path/to/software/simonhmartin/genomics_general/VCF_processing/parseVCF.py -i mydata.vcf.gz --skipIndels | /opt/software/htslib-1.14/bgzip > mydata.geno.gz &



#### ABBA-BABA ANALYSES:
# Prepare a population file: tab-separated, two columns (sample name and population name), no headers

# Define on the command line your P1, P2, P3 (introgressor) and the outgroup
# Define window size in bp (-w), stepsize (-s), and minimum number of SNPs per window (-m)
python /path/to/software/simonhmartin/genomics_general/VCF_processing/ABBABABAwindows.py -g mydata.geno.gz -f phased -o mydata.output.csv -w 100000 -m 100 -s 100000 -P1 EastEurope1 -P2 EastEurope2 -P3 Yarkand -O Altai -T 10 --minData 0.5 --popsFile Brown135_popfile2.txt --writeFailedWindows &


# output file contains the following columns:
# scaffold
# windowstart
# windowend
# nsites	
# nsites_used:	number of sites used to compute statistics (biallelic SNPs)
# ABBA:		pseudo count of ABBA sites (including polymorphic sites, see Martin et al. 2015 Equation 2)
# BABA:		Pseudo count of BABA sites (including polymorphic sites, see Martin et al. 2015 Equation 3)
# D:		D statistic (see Durand et al. 2011 Equation 2)
# fd:		fd admixture estimation (See Martin et al. 2015 Equation 6)
# fdM:		Malinsky's modified statistic, fdM to accomodate admixture between either P1 and P3 or P2 and P3 (See Malinsky et al. 2015 Supplementart Material Page 8)




#### ONCE RUN IS FINISHED, CLOSE CONDA ENVIRONMENT:
conda deactive conenv3.7




### CALCULATE SIGNIFICANCE IN R:

# calculate genome-wide D and the Z-score in R:
# The z-score estimates whether the D-score is significantly different from zero, by calculating how different the D-score is from zero in units of standard error. 
# Based on one comparison, Z-score produces very similar conclusions to jackknifing.

# x		<- read.table("mydata.output.csv",header=TRUE,sep=",")
# mymean	<- mean(x$D,na.rm=TRUE)
# mysd	<- sd(x$D,na.rm=TRUE)
# mystderr	<- mysd/sqrt(nrow(x[is.finite(x$D),]))
# z			<- mymean/mystderr

# plot:
# hist(x$D,breaks=50)
# abline(v=mymean,col="red") 

 

