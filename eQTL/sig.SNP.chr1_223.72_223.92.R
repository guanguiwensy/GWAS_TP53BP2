gwas_0.001_chr1_223.72_223.92 <- results_as[which(results_as$CHR == 1 & results_as$P <= 0.01 & results_as$BP >= 223780000 & results_as$BP <= 223920000),]

pdf("sig.SNP.chr1_223.72_223.92.pdf",width = 15, height = 2)
plot(gwas_0.001_chr1_223.72_223.92$BP,gwas_0.001_chr1_223.72_223.92$CHR,pch="|",xlim = c(223780000,223920000))
dev.off()