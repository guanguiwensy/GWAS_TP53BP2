library("qqman")

results_as <- read.table("assoc_results.assoc", head=TRUE)

color_set <- c("#6675B5","#E78AC3","#A6D854","#FFD92F","#E5C494","#66C2A5","#FC8D62")

par(cex=0.8)

results_as_new <- results_as[results_as$P < 0.1,]
rs7519753_2000000_sig <- results_as_new[which(results_as_new$CHR == 1 & 
                                            results_as_new$BP < 224864224 & 
                                            results_as_new$BP > 222864224 & 
                                            results_as_new$SNP!="."),]

rs7519753_2000000_sig <- c(rs7519753_2000000_sig$SNP)
#Manhattan plot
pdf("manhattan.pdf",height=6,width = 11)
manhattan(results_as_new, col = color_set,suggestiveline = T,
          genomewideline=T,highlight = rs7519753_2000000_sig)
dev.off()

#QQ plot
png("qq.png")
qq(results_as$P, main = "Q-Q plot of GWAS p-values : log")
dev.off()

#inflation factor
p_value <- results_as$P

z <- qnorm(p_value/2)

lambda <- round(median(z^2,na.rm=T)/0.454,3)

lambda