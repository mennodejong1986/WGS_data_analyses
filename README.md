# Collection  of scripts for processing of WGS data

1. Quality control of fastq-files: FASTQ_runfastp.sh
2. Map short reads to reference genome: FASTQ2BAM_mapping_highway.sh
3. Optionally extract reads which mapped to x-chrom and y-chrom: BAM_selectionregion.sh
4. Genotype calling: BAM2VCF_run_mpileup_parallel_highway.sh
5. Obtain pairwise distances for all individuals in vcf-file: VCF_calcdist.sh
6. Obtain sliding-window heterozygosity-estimates for all individuals in vcf-file: VCF_darwindow.sh
7. Population-genetic analyses with R package SambaR, using output of step 4 and 5: VCF_PopGen_SambaR.txt 
