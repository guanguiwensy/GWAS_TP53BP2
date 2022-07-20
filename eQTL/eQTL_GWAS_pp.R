eQTL_GWAS <- read.table("Liver.allpairs.txt.hg19.100000000.txt.hg38.prepaired.txt.0.001.txt.head.txt",header = T,sep = " ")

eQTL_GWAS$gene <- substr(eQTL_GWAS$gene_id,1,15)

protein_coding <- read.table("protein_coding.csv",header = T, sep = "\t")

eQTL_GWAS_protein_coding <- merge(eQTL_GWAS,protein_coding, by.x = "gene", by.y = "Gene.stable.ID")

eQTL_GWAS_protein_coding <- eQTL_GWAS_protein_coding[order(eQTL_GWAS_protein_coding$pval_nominal),]

gwas_0.001_chr1_223.72_223.92 <- results_as[which(results_as$CHR == 1 & results_as$P <= 0.01 & results_as$BP >= 223780000 & results_as$BP <= 223920000),]

pdf("sig.SNP.chr1_223.72_223.92.pdf",width = 15, height = 2)
plot(gwas_0.001_chr1_223.72_223.92$BP,gwas_0.001_chr1_223.72_223.92$CHR,pch="|",xlim = c(223780000,223920000))
dev.off()

gwas_0.001_chr1_223.72_223.92_enhance <- results_as[which(results_as$CHR == 1 & results_as$P <= 0.01 & results_as$BP >= 223886000 & results_as$BP <= 223889000),]


