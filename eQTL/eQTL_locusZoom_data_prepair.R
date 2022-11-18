library(dplyr)

results_as <- read.table("assoc_results.assoc.fisher.txt", head=TRUE)

results_as_new <- results_as[results_as$P < 0.01,]
rs7519753_2000000_sig <- results_as_new[which(results_as_new$CHR == 1 & 
                                                  results_as_new$BP < 224864224 & 
                                                  results_as_new$BP > 222864224 & 
                                                  results_as_new$SNP!="."),]

r2 <- read.csv("rs7519753_2000000_sig.LD.r2.csv")

r2 <- r2[r2$r2 >0.6,]

rs7519753_2000000_sig <- rs7519753_2000000_sig[rs7519753_2000000_sig$SNP %in% r2$RS_number,]

rs7519753_2000000_eQTL <- read.table("Liver.allpairs.txt.hg38.prepaired.rs7519753_2000000.txt",header=T,sep = "")

rs7519753_2000000_eQTL <- rs7519753_2000000_eQTL |> mutate(CHR = t(data.frame(strsplit(chr.pos,split = ":")))[,1],
                                                            BP = t(data.frame(strsplit(chr.pos,split = ":")))[,2],
                                                            A1 = t(data.frame(strsplit(variant_id,split = "_")))[,4],
                                                            A2 = t(data.frame(strsplit(variant_id,split = "_")))[,3],
                                                            P  = pval_nominal,
                                                            SLOPE = slope,
                                                            )

rs7519753_2000000_eQTL <- rs7519753_2000000_eQTL[which(rs7519753_2000000_eQTL$gene_id == "ENSG00000143514.12"),]

rs7519753_2000000_eQTL_plink_format <- rs7519753_2000000_eQTL[,c(11:16)]

rs7519753_2000000_eQTL_plink_format <- rs7519753_2000000_eQTL_plink_format[rs7519753_2000000_eQTL_plink_format$BP %in% rs7519753_2000000_sig$BP,]

write.table(rs7519753_2000000_eQTL_plink_format,"rs7519753_2000000_eQTL_plink_format.txt",row.names = F,col.names = T,sep = "\t",quote = F)
write.table(rs7519753_2000000_sig,"rs7519753_2000000_sig.r2.csv",row.names = F,col.names = T,sep = "\t",quote = F)
