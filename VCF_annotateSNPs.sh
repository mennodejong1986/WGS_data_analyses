#!/bin/bash
# This script is to annotate SNPs using the software SnpEff

####################################
SNPEFF=/home/mdejong/software/snpEff
PREFIX=Brownbear

MYVCF=mydata.filtered.thinned.vcf
REFERENCE=genome.fa
GFF=myannotation.gff
GENES=mygenes.nucl.fa
PROTEINS=mygenes.aa.fa
####################################


# prepare files:
echo "Copying files..."
mkdir ${SNPEFF}/data/${PREFIX}
# mkdir ${SNPEFF}/data/genomes
cp $REFERENCE ${SNPEFF}/data/genomes/${PREFIX}.fa
cp $GFF ${SNPEFF}/data/${PREFIX}/genes.gff
cp $GENES ${SNPEFF}/data/${PREFIX}/sequences.fa
cp $PROTEINS ${SNPEFF}/data/${PREFIX}/protein.fa
echo "${PREFIX}.genome : ${PREFIX}" >> ${SNPEFF}/snpEff.config

# build database:
echo "Building database..."
java -jar ${SNPEFF}/snpEff.jar build -gff3 -v ${PREFIX} &

### SnpEff annotates variant sites only. If we want to annotate all sites (including monomorphic sites), we have to change fifth column from . to A,C,G,T (see awk command below)
# By listing as alternative allele all possible alternate alleles (A,C,G,T), SnpEff will list all possible affects.  
zgrep -B1000000 -a -m 1 '#CHROM' allsites.vcf.gz > input.header.txt
zcat allsites.vcf.gz |  grep -v '#' | awk -F $'\t' '{OFS = FS} { if ($5 == ".") $5="A,C,G,T"; print $0 }' > input.data.txt
cat input.header.txt input.data.txt | gzip > allsites.markedpoly.vcf.gz
rm input.header.txt input.data2.txt
# Note, if necessary, afterwards 'A,C,G,T' can be replaced with . again.

# annotate snps:
echo "Annotating snps..."
java -Xmx4g -jar ${SNPEFF}/snpEff.jar -c ${SNPEFF}/snpEff.config -v $PREFIX mydata.filtered.thinned.vcf > myannotated.vcf &
#java -Xmx4g -jar /home/mdejong/software/snpEff/snpEff.jar -c /home/mdejong/software/snpEff/snpEff.config -v Brownbear allsites.markedpoly.vcf.gz | gzip > annotated.allsites.vcf.gz &


# extracting info about functional snps, to be imported by SambaR:
wait
grep -v '#' myannotated.vcf > mydata.txt
cut -f 2 -d '|' mydata.txt > myeffects.txt
paste myeffects.txt mydata.txt > mydata2.txt

awk '$1 != "downstream_gene_variant" && $1 != "upstream_gene_variant" && $1 != "intergenic_region" && $1 != "intron_variant" && $1 != "splice_region_variant&intron_variant" && $1 != "splice_donor_variant&intron_variant"' mydata2.txt > myfunctional.txt &
sed 's/|,/%/g' myfunctional.txt | cut -f1 -d '%' | sed 's/;ANN/%ANN/g' | cut -f2 -d '%' | cut -f 1 | sed 's/||*/|/g' > myfunctional2.txt
cut -f 2,3,4,9,10,11,13,14 -d '|' myfunctional2.txt | sed 's/|/\t/g' > myfunctional3.txt
 
cut -f 2,3 myfunctional.txt > mysnppositions.txt
paste mysnppositions.txt myfunctional3.txt > snp_ann.txt
sed -i '1 i\chrom\tpos\tclass\timpact\tgene\texon\tnucl_mut\taa_mut\tcodon_nr\taa_nr' snp_ann.txt
rm mysnppositions.txt mydata.txt mydata2.txt myeffects.txt myfunctional3.txt myfunctional2.txt myfunctional3.txt

echo "Done :-)"

# synonymous sites monomorphic:
# zgrep 'synonymous' annotated.allsites.vcf.gz | grep -v 'HIGH' | grep -v 'MODERATE' | head -10
zgrep -v '#' annotated.allsites.vcf.gz | head -n 1000000 | grep 'synonymous' | awk '$5=="A,C,G,T"' | grep -v 'HIGH' | grep -v 'MODERATE' | wc -l
567
# synonymous sites polymorphic:
zgrep -v '#' annotated.allsites.vcf.gz | head -n 1000000 | grep 'synonymous' | awk '$5!="A,C,G,T"' | wc -l
35

# low impact monomorphic:
zgrep -v '#' annotated.allsites.vcf.gz | head -n 10000000 | grep 'LOW' | awk '$5=="A,C,G,T"' | grep -v 'HIGH' | grep -v 'MODERATE' | wc -l
18067
# low impact polymorphic:
zgrep -v '#' annotated.allsites.vcf.gz | head -n 10000000 | grep 'LOW' | awk '$5!="A,C,G,T"' | grep -v 'HIGH' | grep -v 'MODERATE' | wc -l
940

# high impact monomorphic:
zgrep -v '#' annotated.allsites.vcf.gz | head -n 10000000 | grep 'HIGH' | awk '$5=="A,C,G,T"' | grep -v 'LOW' | grep -v 'MODERATE' | wc -l
1886
# high impact polymorphic:
zgrep -v '#' annotated.allsites.vcf.gz | head -n 10000000 | grep 'HIGH' | awk '$5!="A,C,G,T"' | grep -v 'LOW' | grep -v 'MODERATE' | wc -l
5

# dN/dS:
(5/1886)/(35/567)


# classes:
# downstream_gene_variant
# intergenic_region
# intron_variant
# missense_variant
# missense_variant&splice_region_variant
# splice_donor_variant&intron_variant
# splice_region_variant&intron_variant
# splice_region_variant&synonymous_variant
# start_lost
# stop_gained
# synonymous_variant
# upstream_gene_variant

