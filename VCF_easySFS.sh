# See also:
# https://github.com/isaacovercast/easySFS

# VCF to SFS with easySFS:
conda create -n easySFS & conda activate easySFS
conda install -c bioconda dadi pandas
git clone https://github.com/isaacovercast/easySFS.git
cd easySFS
chmod 777 easySFS.py

# Needed is a populations file:
# sample1 pop1
# sample2 pop1
# sample3 pop2
# sample4 pop2
# etc.

# Converting VCF to SFS is a 2 step process. 
# The first step (preview) is to run a preview to identify the values for projecting down each population. In other words: reduce the sample size to a number for which there is no missing data.   
# The second step (convert) is to actually do the conversion specifying the projection values. 


# STEP 1. PREVIEW:
./easySFS.py -i input.vcf -p pops_file.txt -a --preview --total-length 10000
# see below for explanation of flags.

# This will output something like:
# pop1
# (2, 45.0)   (3, 59.0)   (4, 58.0)   (5, 49.0)     
# pop2
# (2, 68.0)   (3, 96.0)   (4, 106.0)  (5, 110.0)  

# Each column is the number of (haploid!) samples in the projection and the number of segregating sites at that projection value. 
# The dadi manual recommends maximizing the number of segregating sites.
# However, if you have lots of missing data then you might have to balance number of segregating sites against number of samples, in order to avoid downsampling too far.



# STEP 2. CONVERT:
./easySFS.py -i input.vcf -p pops_file.txt -a --proj=4,5 --total-length 10000
# --window-bp:	Select SNPs based on window size in base pairs (optional, to trim down)
# --total-length:	total sequence length of input data (for accurate zero bin)
# -a: 		include all sites in the vcf file
# --proj:	list of values for projecting populations down to different sample sizes. 
#			In this example for 2 populations, the first population will be estimates for 4 (haploid!!) individuals, and the second population for 5 individuals.
#			The preview tells us that this means: calculate SFS over 58 sites for population 1, and to calculate SFS over 106 sites for population 2. 



#### HOW TO USE OUTPUT FOR STAIRWAY PLOT?

# How to interpret the output?
# Say that you have a dataset of 10 individuals, projected down to 6 (haploid!) individuals.
# The folded sfs-vector outputted by easysfs will look something like this (note that there will be two rows):

# 1000 100 95 90 0 0 0 
# 1 0 0 0 1 1 1

# IMPORTANT!! 
# The second row, with the zeros and ones, is not part of the sfs-vector.
# Instead, it is another vector, a boolean vector, which indicates which bins to select.
# In the words of the author: 
# 'The zeros and ones are a bitmask to indicate which bins should be considered as data.' 
# 'In an unfolded sfs this whole row of output will be 1s, in a folded sfs the row will contain 0s and 1s in roughly equal proportion, as a result of the folding procedure.'
# 0 means: select. 
# 1 means: do not select.
# So in the above example, the final sfs-vector (to be selected and used as input for Stairwayplot), becomes: 100 95 90


