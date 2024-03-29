#!/bin/bash

#$anno means dbSNP annotate vcf file;
#SNP_quality_control

bcftools view --types snps -m 2 -M 2 --exclude-uncalled \
-O z -o merge.filted.vcf.gz  merge.vcf.gz

bcftools annotate -a $anno -c ID -O z -o merge.filted.annotate.vcf.gz  merge.filted.vcf.gz