# Collection of scripts for processing of WGS data
The file 'GenotypeCalling.routemap.pdf' contains a schematic overview of the overall workflow.

1. Quality control of fastq-files: FASTQ_runfastp.sh
2. Map short reads to reference genome: FASTQ2BAM_bwa_HIGHWAY.sh
3. Optionally extract reads which mapped to x-chrom and y-chrom: BAM_selectregion.sh
4. Genotype calling: BAM2VCF_bcftools_HIGHWAY.sh
5. Obtain pairwise distances for all individuals in vcf-file: VCF_calcdist.sh
6. Obtain sliding-window heterozygosity-estimates for all individuals in vcf-file: VCF_darwindow.sh
7. Population-genetic analyses with R package SambaR, using output of step 4 and 5: VCF_PopGen_SambaR.txt
8. For other popgen-analyses, use remaining scripts starting with the prefix 'VCF'.

# Default filter settings
fastq-files:  
minimum base quality: 15        fastp (FASTQ_runfastp.sh)  
  
bam-files:  
mimimum mapping quality: 20     samtools view -q 20 (FASTQ2BAM_bwa_highway.sh)  
minimum alignment score: 100    samtools view -e '[AS]>=100' (FASTQ2BAM_bwa_highway.sh)  
properly pairs reads only       samtools view -f 0x2 (FASTQ2BAM_bwa_highway.sh)  
  
vcf-files:  
min. depth per genotype: 5      bcftools filter --set-GTs . -e "FMT/DP<5" (BAM2VCF_bcftools_HIGHWAY.sh)  
min. depth per site:            customised value based on observed depth distribution (see BAM2VCF_plotdepth_inR.txt)  
max. depth per site:            idem  
no indels                       bcftools view --exclude-types indels  
  
Because filters may introduce bias, the pipeline does not contain filters on site quality, minor allele frequency, missingness per site, or linkage disequilibrium.  

