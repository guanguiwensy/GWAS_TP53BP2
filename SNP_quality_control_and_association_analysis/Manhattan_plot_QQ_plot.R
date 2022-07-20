library("qqman")

results_as <- read.table("assoc_results.assoc", head=TRUE)

color_set <- c("#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#66C2A5","#FC8D62")

par(cex=0.8)


#Manhattan plot
pdf("manhattan.pdf")
manhattan(results_as, col = color_set,annotatePval = 0.01,suggestiveline = F,genomewideline=F)
dev.off()

#QQ plot
pdf("manhattan.pdf")
qq(results_as$P, main = "Q-Q plot of GWAS p-values : log")
dev.off()

#inflation factor
p_value <- results_as$P

z <- qnorm(p_value/2)

lambda <- round(median(z^2,na.rm=T)/0.454,3)

lambda