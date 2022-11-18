eQTL_GWAS <- read.table("Liver.allpairs.txt.hg38.prepaired.txt.0.0001.txt",header = F,sep = " ")

colnames(eQTL_GWAS) <- c("chr:pos",
                         "gene_id",
                         "variant_id",
                         "tss_distance",
                         "ma_samples",
                         "ma_count",
                         "maf",
                         "pval_nominal",
                         "slope",
                         "slope_se")

eQTL_GWAS$gene <- substr(eQTL_GWAS$gene_id,1,15)

protein_coding <- read.table("protein_coding.csv",header = T, sep = "\t")

eQTL_GWAS_protein_coding <- merge(eQTL_GWAS,protein_coding, by.x = "gene", by.y = "Gene.stable.ID")

write.table(eQTL_GWAS_protein_coding,"eQTL_GWAS_protein_coding.txt",sep = "\t", col.names = T, row.names = F, quote = F)