results_as <- read.table("./gwas/assoc_results.assoc", head=TRUE)

rs7519753_2000000 <- results_as[which(results_as$CHR == 1 & results_as$BP < 224864224 & results_as$BP > 222864224),]

rs7519753_2000000_pos <- paste("chr",rs7519753_2000000[,1],":",rs7519753_2000000[,3],sep="")

write.table(rs7519753_2000000_pos,"rs7519753_2000000_pos.txt",row.names = F,col.names = F,quote = F)

rs7519753_2000000_eQTL <- read.table("Liver.allpairs.txt.hg38.prepaired.txt.0.0001.txt",header=T,sep = "")

rs7519753_2000000_eQTL <- rs7519753_2000000_eQTL %>% mutate(CHR = t(data.frame(strsplit(chr.pos,split = ":")))[,1],
                                                            BP = t(data.frame(strsplit(chr.pos,split = ":")))[,2],
                                                            A1 = t(data.frame(strsplit(variant_id,split = "_")))[,4],
                                                            A2 = t(data.frame(strsplit(variant_id,split = "_")))[,3],
                                                            P  = pval_nominal,
                                                            SLOPE = slope,
                                                            )

rs7519753_2000000_eQTL <- rs7519753_2000000_eQTL[which(rs7519753_2000000_eQTL$gene_id == "ENSG00000143514.12"),]

rs7519753_2000000_eQTL_plink_format <- rs7519753_2000000_eQTL[,c(11:16)]

write.table(rs7519753_2000000_eQTL_plink_format,"rs7519753_2000000_eQTL_plink_format.txt",row.names = F,col.names = T,sep = " ",quote = F)