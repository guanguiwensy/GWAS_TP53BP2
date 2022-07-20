

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

gunzip hg19ToHg38.over.chain.gz

awk -v FS=" " -v OFS="," 'NR>1 {$1=$1;print}' $id >$id.all.txt
awk 'NR>1 {print $2}' $id |awk -v FS="_" -v OFS=" " '{$1=$1;print "chr"$1" "$2" "$2}' - > $id.bed.txt
paste $id.bed.txt $id.all.txt >$id.hg19.txt



sed -n '1,20000000p' Liver.allpairs.txt.hg19.txt >Liver.allpairs.txt.hg19.20000000.txt


tem
arr=($tem)
for i in $(seq 0 6); do ./liftOver ${arr[$i]} hg19ToHg38.over.chain ${arr[$i]}.hg38.txt ${arr[$i]}.unmap & done



for i in $(seq 0 6); do awk '{print $1":"$2","$4}' ${arr[$i]}.hg38.txt|awk -v FS="," -v OFS=" " 'BEGIN{print "chr:pos","gene_id","variant_id","tss_distance","ma_samples","ma_count","maf","pval_nominal","slope","slope_se"}{$1=$1;print}' - >${arr[$i]}.hg38.prepaired.txt & done

awk '{print $11}' results_as_0.001.txt >id_0.01.txt

for i in $(seq 0 6); do perl id_extract_gtex.pl id_0.001.txt ${arr[$i]}.hg38.prepaired.txt >${arr[$i]}.hg38.prepaired.txt.0.0001.txt & done

for i in $(seq 1 6); do cat ${arr[$i]}.hg38.prepaired.txt.0.001.txt >> ${arr[0]}.hg38.prepaired.txt.0.001.txt ;done

head -n 1 Liver.allpairs.txt.hg19.20000000.txt.hg38.prepaired.txt > ${arr[0]}.hg38.prepaired.txt.0.001.txt.head.txt

cat ${arr[0]}.hg38.prepaired.txt.0.001.txt >>${arr[0]}.hg38.prepaired.txt.0.001.txt.head.txt

Rscript eQTL_protein_gene.R


