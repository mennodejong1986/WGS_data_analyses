#!/bin/bash
# convert haploid vcf to diploid vcf:

###################
#MYVCF=Brown135_xchrom.mysnps.thinned.2000.vcf.gz
MYVCF=NorthAmerica112_xchrom.mysnps.thinned.2000.vcf.gz
###################


echo "Converting haploid data to diploid data..."
zcat $MYVCF | grep '#' > mydiploid.vcf.header
zcat $MYVCF | grep -v '#' | sed 's|0:|0/0:|g' | sed 's|1:|1/1:|g' | sed 's|\.:|\./\.:|g' > mydiploid.vcf.data

# in case input data is haplodiploid (such as X-chrom data), we need to make a correction:
echo "Making corrections in case of haplodiploid data..."
grep -v '#' mydiploid.vcf.data | sed 's|0/0/0:|0/0:|g' | sed 's|0/1/1:|0/1:|g' | sed 's|1/1/1:|1/1:|g' | sed 's|\./\./\.:|\./\.:|g' > mydiploid.vcf.tmp.data
mv mydiploid.vcf.tmp.data mydiploid.vcf.data

echo "Combining..."
cat mydiploid.vcf.header mydiploid.vcf.data | gzip > mydiploid.vcf.gz
rm mydiploid.vcf.header mydiploid.vcf.data

echo "Created a file called mydiploid.vcf.gz."
