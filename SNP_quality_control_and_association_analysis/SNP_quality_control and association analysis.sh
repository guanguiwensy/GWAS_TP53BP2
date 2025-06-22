#data prepair
bgzip -@ 40 merge.vcf  && bcftools index -t --threads 40 merge.vcf.gz

bcftools view --types snps -m 2 -M 2 -q 0.01:minor --threads 40 --exclude-uncalled\
- O z -o merge.vcf.filted.gz merge.vcf.gz

bcftools annotate -a /home/guanguiwen/data1/reference/human/hg38/dbsnp_146.hg38.vcf.gz \
-c ID -O z --threads 40 -o merge.vcf.filted.annotate.gz merge.vcf.filted.gz

plink --vcf merge.vcf.filted.annotate.gz  --make-bed --out 100_HBsAg_GWAS


#Filter SNPs that do not meet the requirements
plink --bfile 100_HBsAg_GWAS --geno 0.2 --make-bed --out 100_HBsAg_GWAS_2

plink --bfile 100_HBsAg_GWAS_2 --mind 0.2 --make-bed --out 100_HBsAg_GWAS_3

plink --bfile 100_HBsAg_GWAS_3 --geno 0.05 --make-bed --out 100_HBsAg_GWAS_4

plink --bfile 100_HBsAg_GWAS_4 --mind 0.02 --make-bed --out 100_HBsAg_GWAS_5

#Gender check
plink --bfile 100_HBsAg_GWAS_5 --check-sex

plink --bfile 100_HBsAg_GWAS_5 --chr 1-22 --make-bed --out 100_HBsAg_GWAS_7

#MAF check and filter SNPs that MAF < 0.05
plink --bfile 100_HBsAg_GWAS_7 --freq --out MAF_check

Rscript --no-save MAF_check.R

#Filter SNPs that MAF < 0.05
plink --bfile 100_HBsAg_GWAS_7 --maf 0.05 --make-bed --out 100_HBsAg_GWAS_8

#HWE check and filter SNPs that HWE < 1e-10

plink --bfile 100_HBsAg_GWAS_8 --hardy

awk '{ if ($9 <0.00001) print $0 }' plink.hwe>plinkzoomhwe.hwe

Rscript --no-save hwe.R

plink --bfile 100_HBsAg_GWAS_8 --hwe 1e-10 --make-bed --out 100_HBsAg_GWAS_9

plink --bfile 100_HBsAg_GWAS_9 --exclude inversion.txt --range --indep-pairwise 50 5 0.2 --out indepSNP

plink --bfile 100_HBsAg_GWAS_9 --extract indepSNP.prune.in --het --out R_check

Rscript --no-save check_heterozygosity_rate.R

Rscript --no-save heterozygosity_outliers_list.R

sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt

plink --bfile 100_HBsAg_GWAS_9 --remove het_fail_ind.txt --make-bed --out 100_HBsAg_GWAS_10

#For each pair of 'related' individuals with a pihat > 0.2, we remove the individual with the lowest call rate
plink --bfile 100_HBsAg_GWAS_10 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2

plink --bfile 100_HBsAg_GWAS_10 --filter-founders --make-bed --out 100_HBsAg_GWAS_11

plink --bfile 100_HBsAg_GWAS_11 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2_in_founders

plink --bfile 100_HBsAg_GWAS_11 --missing


echo "1-55 1" >0.2_low_call_rate_pihat.txt

plink --bfile 100_HBsAg_GWAS_11 --remove 0.2_low_call_rate_pihat.txt --make-bed --out 100_HBsAg_GWAS_12



#The 1000 genome phase3 data prepairation.

plink2 --bfile 1kg.all_hg38 --geno 0.2 --allow-no-sex --allow-extra-chr --chr 1-22 --make-bed --out 1kG_MDS

plink --bfile 1kG_MDS --mind 0.2 --allow-no-sex --make-bed --out 1kG_MDS2

plink --bfile 1kG_MDS2 --geno 0.02 --allow-no-sex --make-bed --out 1kG_MDS3

plink --bfile 1kG_MDS3 --mind 0.02 --allow-no-sex --make-bed --out 1kG_MDS4

plink --bfile 1kG_MDS4 --maf 0.05 --allow-no-sex --make-bed --out 1kG_MDS5

awk '{print$2}' 100_HBsAg_GWAS_12.bim > 100_HBsAg_GWAS_12_SNPs.txt

plink2 -bfile 1kG_MDS5 --recode vcf-iid --output-chr chrM --allow-extra-chr --chr 1-22 --keep sample_keep -out 1kG_MDS5

bgzip -@ 40 1kG_MDS5.vcf  && bcftools index -t --threads 40 1kG_MDS5.vcf.gz

#Annotation of 1000 genome SNPs using dbSNP
bcftools annotate -a /home/guanguiwen/data1/reference/human/hg38/dbsnp_146.hg38.vcf.gz \
-c ID -O z --threads 40 -o 1kG_MDS5.vcf.annotate.gz 1kG_MDS5.vcf.gz

plink --vcf 1kG_MDS5.vcf.annotate.gz  --make-bed --out 1kG_MDS5.1

mv 1kG_MDS5.1.bim 1kG_MDS5.bim

plink --bfile 1kG_MDS5 --extract 100_HBsAg_GWAS_12_SNPs.txt --make-bed --out 1kG_MDS6

#SNPs detected in both 1000 genomes and our own samples were retained, and the SNP types were guaranteed to be consistent.
awk '{print$2}' 1kG_MDS6.bim > 1kG_MDS6_SNPs.txt


plink --bfile 100_HBsAg_GWAS_12 --extract 1kG_MDS6_SNPs.txt --recode --make-bed --out 100_HBsAg_GWAS_12_MDS


awk '{print$2,$4}'  100_HBsAg_GWAS_12_MDS.map > build_100_HBsAg_GWAS_map.txt

plink --bfile 1kG_MDS6 --update-map build_100_HBsAg_GWAS_map.txt --make-bed --out 1kG_MDS7

awk '{print$2,$5}' 1kG_MDS7.bim > 1kg_ref-list.txt

plink --bfile 100_HBsAg_GWAS_12_MDS --reference-allele 1kg_ref-list.txt --make-bed --out 100_HBsAg_GWAS-adj

awk '{print$2,$5,$6}' 1kG_MDS7.bim > 1kGMDS7_tmp

awk '{print$2,$5,$6}' 100_HBsAg_GWAS-adj.bim > 100_HBsAg_GWAS-adj_tmp

sort 1kGMDS7_tmp 100_HBsAg_GWAS-adj_tmp |uniq -u > all_differences.txt

awk '{print$1}' all_differences.txt | sort -u > flip_list.txt

plink --bfile 100_HBsAg_GWAS-adj --flip flip_list.txt --reference-allele 1kg_ref-list.txt --make-bed --out corrected_100_HBsAg_GWAS

awk '{print$2,$5,$6}' corrected_100_HBsAg_GWAS.bim > corrected_100_HBsAg_GWAS_tmp

sort 1kGMDS7_tmp corrected_100_HBsAg_GWAS_tmp |uniq -u  > uncorresponding_SNPs.txt

awk '{print$1}' uncorresponding_SNPs.txt | sort -u > SNPs_for_exlusion.txt

plink --bfile corrected_100_HBsAg_GWAS --exclude SNPs_for_exlusion.txt --make-bed --out 100_HBsAg_GWAS_12_MDS2

plink --bfile 1kG_MDS7 --exclude SNPs_for_exlusion.txt --make-bed --out 1kG_MDS8


#Merge our sample and 1000 genome
plink --bfile 100_HBsAg_GWAS_12_MDS2 --bmerge 1kG_MDS8.bed 1kG_MDS8.bim 1kG_MDS8.fam --allow-no-sex --make-bed --out MDS_merge2

plink --bfile MDS_merge2 --extract indepSNP.prune.in --genome --out MDS_merge2

plink --bfile MDS_merge2 --read-genome MDS_merge2.genome --cluster --mds-plot 10 --out MDS_merge2

awk '{print $1,$2,"OWN"}' 100_HBsAg_GWAS_12_MDS.fam > racefile_own.txt


#MDS plot.
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/20100804.ALL.panel
awk '{print$1,$1,$2}' 20100804.ALL.panel > race_1kG.txt
sed 's/JPT/ASN/g' race_1kG.txt>race_1kG2.txt
sed 's/ASW/AFR/g' race_1kG2.txt>race_1kG3.txt
sed 's/CEU/EUR/g' race_1kG3.txt>race_1kG4.txt
sed 's/CHB/ASN/g' race_1kG4.txt>race_1kG5.txt
sed 's/CHD/ASN/g' race_1kG5.txt>race_1kG6.txt
sed 's/YRI/AFR/g' race_1kG6.txt>race_1kG7.txt
sed 's/LWK/AFR/g' race_1kG7.txt>race_1kG8.txt
sed 's/TSI/EUR/g' race_1kG8.txt>race_1kG9.txt
sed 's/MXL/AMR/g' race_1kG9.txt>race_1kG10.txt
sed 's/GBR/EUR/g' race_1kG10.txt>race_1kG11.txt
sed 's/FIN/EUR/g' race_1kG11.txt>race_1kG12.txt
sed 's/CHS/ASN/g' race_1kG12.txt>race_1kG13.txt
sed 's/PUR/AMR/g' race_1kG13.txt>race_1kG14.txt


cat race_1kG14.txt racefile_own.txt | sed -e '1i\FID IID race' > racefile.txt

awk '{if($1==0) $1=$2;print}' MDS_merge2.mds >MDS_merge2.mds1 &&mv MDS_merge2.mds1 MDS_merge2.mds

Rscript MDS_merged.R


plink --bfile 100_HBsAg_GWAS_12_MDS2 --extract indepSNP.prune.in --genome --out 100_HBsAg_GWAS_12_MDS2

plink --bfile 100_HBsAg_GWAS_12_MDS2 --read-genome 100_HBsAg_GWAS_12_MDS2.genome --cluster --mds-plot 10 --out own

awk '{print $1,$2,$6}' 100_HBsAg_GWAS_12_MDS.fam > raceown.txt

Rscript MDS_own.R

#Association analysis
plink --bfile 100_HBsAg_GWAS_12 --extract indepSNP.prune.in --genome --out 100_HBsAg_GWAS_12

plink --bfile 100_HBsAg_GWAS_12 --read-genome 100_HBsAg_GWAS_12.genome --cluster --mds-plot 10 --out 100_HBsAg_GWAS_12

awk '{print$1, $2, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' 100_HBsAg_GWAS_12.mds > covar_mds.txt

plink --bfile 100_HBsAg_GWAS_12 --assoc --out assoc_results

plink --bfile 100_HBsAg_GWAS_12 --covar covar_mds.txt --logistic --hide-covar --out logistic_results

awk '!/'NA'/' logistic_results.assoc.logistic > logistic_results.assoc_2.logistic

sort -g -k 9 logistic_results.assoc_2.logistic >> logistic_results.assoc_2.logistic.sort

plink --annotate assoc_results.assoc.fisher ranges=glist-hg38  distance --border 1000
