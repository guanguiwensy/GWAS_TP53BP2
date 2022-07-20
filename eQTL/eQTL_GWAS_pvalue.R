
library(magrittr)

library(dplyr)

library(clusterProfiler)

library(org.Hs.eg.db)

library("ggrepel")

gtex=read.table("Liver.allpairs.txt.hg19.20000000.txt.hg38.prepaired.txt.0.0001.txt.head.txt",sep=" ",header=T)

id_split <- gtex$variant_id %>% strsplit(split="_") %>% data.frame() %>% t() %>% as.data.frame()

gtex$hg38id <- paste(gtex[,1],id_split[,3],id_split[,4],sep = "_")

proteincoding <- read.table("protein_coding.csv",sep="\t",header=T)

gtex$gene_id_1 <-  substr(gtex[,2],1,15)

gtex <- merge(gtex,proteincoding,by.x = "gene_id_1",by.y = "Gene.stable.ID")




gwas_0.0001 <- read.table("results_as_0.0001.txt",header = T, sep = " ")

gwas_0.0001$hg38id <- paste("chr",gwas_0.0001[,1],":",gwas_0.0001[,3],"_",gwas_0.0001[,7],"_",gwas_0.0001[,4],sep="")

gwas_eQTL_merge <- merge(gwas_0.0001,gtex,by.x = "hg38id",by.y = "hg38id")

gwas_eQTL_merge <- gwas_eQTL_merge[order(gwas_eQTL_merge$pval_nominal),]

TP53BP2 <- gwas_eQTL_merge[which(gwas_eQTL_merge$Gene.name=="TP53BP2"),]



gwas_eQTL_merge <- gwas_eQTL_merge %>% mutate(abs_clour=case_when(slope > 0 ~ "red",TRUE ~ "blue"))

gwas_eQTL_merge <- gwas_eQTL_merge %>% mutate(highlight=case_when(abs(slope) > 0.15 ~ paste(SNP,":",Gene.name,sep=""),TRUE ~ ""))

plot <- ggplot(gwas_eQTL_merge, aes(x = -log10(P), y = -log10(pval_nominal))) +
               geom_point(data = gwas_eQTL_merge[which(gwas_eQTL_merge$Gene.name != "TP53BP2"),],
                          aes(color = I("skyblue"),
                              size = abs(slope), 
                              alpha = 0.9))+
               geom_point(data = gwas_eQTL_merge[which(gwas_eQTL_merge$Gene.name == "TP53BP2"),],
                          aes(x = -log10(P),
                              y = -log10(pval_nominal),
                              color=I("red"),
                              size=abs(slope),
                              alpha = 0.9))+
               geom_text_repel(aes(label = highlight),
                               size = 3,
                               min.segment.length = 0, 
                               seed = 42, 
                               box.padding = 0.5,
                               max.overlaps = 9,
                               arrow = arrow(length = unit(0.010, "npc")),
                               nudge_x = .15,
                               nudge_y = .5,
                               color = "grey50"
               )



results_as <- read.table("./gwas/assoc_results.assoc", head=TRUE)

rs7519753_2000000 <- results_as[which(results_as$CHR == 1 & results_as$BP < 224864224 & results_as$BP > 222864224),]

rs7519753_2000000_pos <- paste("chr",rs7519753_2000000[,1],":",rs7519753_2000000[,3],sep="")

write.table(rs7519753_2000000_pos,"rs7519753_2000000_pos.txt",row.names = F,col.names = F,quote = F)

rs7519753_2000000_eQTL <- read.table("Liver.allpairs.txt.hg19.20000000.txt.hg38.prepaired.txt.rs7519753.head.txt",header=T,sep = "")

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

eQTL_pos_TP53BP2=ggplot(gwas_eQTL_merge, aes(x = BP, y = -log10(pval_nominal))) +
  geom_point(data = gwas_eQTL_merge[which(gwas_eQTL_merge$Gene.name == "TP53BP2"),],
             aes(color = I("skyblue"),
                 size = abs(slope), 
                 alpha = 0.9))+xlim(223800000,223900000)

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
