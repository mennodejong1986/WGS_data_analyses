# Collection  of scripts for processing of WGS data
The file 'GenotypeCalling.routemap.pdf' contains a schematic overview of the overall workflow.

1. Quality control of fastq-files: FASTQ_runfastp.sh
2. Map short reads to reference genome: FASTQ2BAM_bwa_HIGHWAY.sh
3. Optionally extract reads which mapped to x-chrom and y-chrom: BAM_selectregion.sh
4. Genotype calling: BAM2VCF_bcftools_HIGHWAY.sh
5. Obtain pairwise distances for all individuals in vcf-file: VCF_calcdist.sh
6. Obtain sliding-window heterozygosity-estimates for all individuals in vcf-file: VCF_darwindow.sh
7. Population-genetic analyses with R package SambaR, using output of step 4 and 5: VCF_PopGen_SambaR.txt
8. For other popgen-analyses, use remaining scripts starting with the prefix 'VCF'.
