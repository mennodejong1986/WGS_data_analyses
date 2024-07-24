#!/bin/bash
# script to convert (gzipped) VCF to PED/MAP format (and binary RAW/BIM) using vcftools

###################
VCFTOOLS=/opt/software/vcftools/vcftools_0.1.17/bin/vcftools
PLINK=/home/mdejong/software/plinkv2020-09-21/plink		

MYVCF=mydiploid.vcf.gz		# input gzipped vcf file
PREFIX=Gobi_xchrom	# prefix of output files
###################

echo "Converting vcf to ped/map..."
${VCFTOOLS} --gzvcf ${MYVCF} --plink --out $PREFIX
# an alternative is to use plink: plink --vcf ${MYVCF} --recode --out ${PREFIX}

echo "Editing first column of map file..."
cut -f2 ${PREFIX}.map | cut -f1 -d ':' > mycontigs.txt && cut -f2,3,4 ${PREFIX}.map > mymap.txt && paste mycontigs.txt mymap.txt > ${PREFIX}.map && rm mycontigs.txt mymap.txt

echo "Converting ped/map to raw/bim..."
${PLINK} --file ${PREFIX} --allow-extra-chr --chr-set 95 --make-bed --recode A --out ${PREFIX}

echo "Finished converting vcf file :-)"



