#!/bin/bash

#Downloading GTEx v7 and liftOver chain file.
wget https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/all_snp_gene_associations/Liver.allpairs.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

gunzip Liver.allpairs.txt.gz
gunzip hg19ToHg38.over.chain.gz

#Data prepair for liftOver
$id=Liver.allpairs.txt
awk -v FS=" " -v OFS="," 'NR>1 {$1=$1;print}' $id >$id.all.txt
awk 'NR>1 {print $2}' $id |awk -v FS="_" -v OFS=" " '{$1=$1;print "chr"$1" "$2" "$2}' - > $id.bed.txt
paste $id.bed.txt $id.all.txt >$id.hg19.txt

liftOver $id.hg19.txt hg19ToHg38.over.chain $id.hg38.txt $id.unmap

#Merging of GTEx v7 data and GWAS data
awk '{print $1":"$2","$4}' $id.hg38.txt|awk -v FS="," -v OFS=" " \
'BEGIN{print "chr:pos","gene_id","variant_id","tss_distance","ma_samples","ma_count","maf","pval_nominal","slope","slope_se"}{$1=$1;print}' - >$id.hg38.prepaired.txt

awk '{print $11}' results_as_0.001.txt >id_0.001.txt

perl id_extract_gtex.pl id_0.001.txt $id.hg38.prepaired.txt >$id.hg38.prepaired.txt.0.0001.txt

#Keeping protein coding genes.
Rscript eQTL_protein_gene.R