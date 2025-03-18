#!/bin/bash
# script for sliding window Dxy, pi, He, Fst and ABBA-BABA analyses (D and f statistics)
# Using software developed by Simon Martin [Martin et al. 2015, MBE]

# Online tutorial:
# http://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/abba-baba-statistics/

# Note: the tutorial mentioned that this software runs in a python 2.7 environment, but this seems to be incorrect and leads to errors (such as the one explained below).
# The software requires python3.7. 
# Do not copy the executables of Simon Martin's software to your working directory, because they may depend on other executables. 
# For example, parseVCF.py depends on executable genomics.py. So keep all together in the git clone directory.

# you can either use python or python3.
# python3 can avoid the error for ABBAwindows.py:
# line 64: print("Sorter received result", resNumber, file=sys.stderr) ^ SyntaxError: invalid syntax



### INSTALLATION

# git clone https://github.com/simonhmartin/genomics_general

# conda activate /home/mdejong/software/anaconda3/envs/conenv2.7
# conda install numpy
# cd /home/mdejong/bearproject/snpfiles_december2020/slidingwindow121/testrun




### CONVERT VCF TO GENO FORMAT

# convert input vcf to geno format (format developed by Simon Martin):
# parseVCF.py for single thread and parseVCFs.py (with an 's') for multiple threads:
# If wanted to include a site quality filter, add after --skipIndels: --minQual 30
python3 /home/mdejong/software/simonhmartin/genomics_general/VCF_processing/parseVCF.py -i mydiploid.vcf.gz --skipIndels | /opt/software/htslib-1.18/bgzip > mysnps.missfilter.geno.gz &

# parseVCFs divides the vcf files into sections (1 per thread), and will run the conversion for all contigs specified in the vcf header section
# Note: when running the command, you might receive this error:
# parseVCF.py on line 195: TypeError: split() takes no keyword arguments
# This is likely caused by using the wrong Python environment (Python 2.7 instead of Python 3.7) 
# If so, edit the parseVCF.py script as follows:
# contigDataDict = dict([x.split("=", 1) for x in re.split('<|>', line)[1].split(",")])



### GENERAL OPTIONS FOR SLIDING WINDOW ANALYSES
# --windType {sites,coordinate,predefined}			Type of windows to make
# -w sites, --windSize sites					Window size in bases
# -s sites, --stepSize sites					Step size for sliding window
# -m sites, --minSites sites					Minumum good sites per window
# --overlap sites       					Overlap for sites sliding window
# -D MAXDIST, --maxDist MAXDIST					Maximum span distance for sites window
# T: number of threads




### SLIDING WINDOW HETEROZYGOSITY, PI, FST AND DXY
# BETWEEN POPULATIONS

# popsFile should have two columns, named 'name' and 'pop'.

# heterozygosity:
python /home/mdejong/software/simonhmartin/genomics_general/popgenWindows.py -w 10000 -m 100 -s 5000 -g trial.geno.gz -o trial.output.csv.gz --windType coordinate -f phased -T 1 --popsFile mypopfile2.txt --analysis indHet &
# Dxy:
python /home/mdejong/software/simonhmartin/genomics_general/popgenWindows.py -w 10000 -m 10 -s 5000 -g trial.geno.gz -o trial.output.csv.gz --writeFailedWindows --windType coordinate -f phased -T 1 -p ABC -p HudsonBay -p Ural --popsFile mypopfile2.txt --analysis popPairDist &
# Fst:
python /home/mdejong/software/simonhmartin/genomics_general/popgenWindows.py -w 10000 -m 10 -s 5000 -g trial.geno.gz -o trial.freq.output.csv.gz --writeFailedWindows --windType coordinate -f phased -T 1 -p ABC -p HudsonBay -p Ural --popsFile mypopfile2.txt --analysis popFreq &

# w: window size (bp)
# m: minimum number of usable sites in window (bp)
# s: step size (bp)
# T: number of threads
# --analysis {popFreq,popDist,popPairDist,indPairDist,indHet} [{popFreq,popDist,popPairDist,indPairDist,indHet} ...]
 

# BETWEEN INDIVIDUALS
# individual pair distance matrix:
python /home/mdejong/software/simonhmartin/genomics_general/distMat.py -g trial.geno.gz -w 10000 -m 100 -s 5000 -f phased --windType coordinate -o trial.output.dist &




### GENOME FRAGMENTS PHYLOGENY

# sliding window trees with phyml and/or raxml (phymlWindows.py and raxmlWindows.py):
python /home/mdejong/software/simonhmartin/genomics_general/phylo/raxml_sliding_windows.py -h
# In python3.7 I received a syntax error at line 329, which as suggested by the error message I managed to solve by replacing "\nDone." with ("\nDone.")
# Then I received an error that the software was not able to located the module 'genomics.py'. 
# This I solved by moving the raxml_sliding_windows.py script from genomics_general/phylo/ to the parent directory genomics_general, which contains the module genomics.py:
mv /home/mdejong/software/simonhmartin/genomics_general/phylo/raxml_sliding_windows.py /home/mdejong/software/simonhmartin/genomics_general/raxml_sliding_windows.py
python /home/mdejong/software/simonhmartin/genomics_general/raxml_sliding_windows.py -T 1 -g mysnps.geno.gz --prefix trial.phymltree.w50 -w 100000 -M 1000 -S 100000 --windType coordinate --model GTR --include myscaffolds.txt --individuals CentralRussia1,Ural2,Rumania1,Rumania2,ABC1,TorontoZoo --outgroup TorontoZoo --raxml /home/mdejong/software/RAXML/standard-RAxML/raxmlHPC &




### ABBA-BABA
# sliding window ABBA-BABA (runs in python3):
# -P1  
# -P2 
# -P3 	introgressor
# -O 	outgroup
# python /home/mdejong/software/simonhmartin/genomics_general/ABBABABAwindows.py -h

# popsFile should have two columns, named 'name' and 'pop'.

python /home/mdejong/software/simonhmartin/genomics_general/ABBABABAwindows.py -g trial.geno.gz -f phased -o trial.Dstats.csv -w 20000 -m 100 -s 10000 -P1 ABC -P2 HudsonBay -P3 polar -O Europe -T 10 --minData 0.5 --popsFile mypopfile2.txt --writeFailedWindows &
# sed 's/,/\t/g' mysnps.Dstats.csv > mysnps.Dstats.txt

python /home/mdejong/software/simonhmartin/genomics_general/ABBABABAwindows.py -g Brown135.mysnps.missfilter.geno.gz -f phased -o Hokkaido_Amur.100000.csv -w 100000 -m 100 -s 100000 -P1 Hokkaido -P2 Amur -P3 polar -O Black -T 10 --minData 0.5 --popsFile Brown135_popfile2.txt --writeFailedWindows &

# calculate D and Z-score (significantly different from zero?) in R:
# x			<- read.table("Hokkaido_Amur.20000.csv",header=TRUE,sep=",")
# mymean	<- mean(x$D,na.rm=TRUE)
# mysd		<- sd(x$D,na.rm=TRUE)
# mystderr	<- mysd/sqrt(nrow(x[is.finite(x$D),]))
# z			<- mymean/mystderr
# Based on one comparison, Z-score produces very similar conclusions to jackknifing. 




#### OUTPUT FILE ####

# output file:
# scaffold
# windowstart
# windowend
# nsites	
# nsites_used:	number of sites used to compute statistics (biallelic SNPs)
# ABBA:			pseudo count of ABBA sites (including polymorphic sites, see Martin et al. 2015 Equation 2)
# BABA:			Pseudo count of BABA sites (including polymorphic sites, see Martin et al. 2015 Equation 3)
# D:			D statistic (see Durand et al. 2011 Equation 2)
# fd:			fd admixture estimation (See Martin et al. 2015 Equation 6)
# fdM:			Malinsky's modified statistic, fdM to accomodate admixture between either P1 and P3 or P2 and P3 (See Malinsky et al. 2015 Supplementart Material Page 8)