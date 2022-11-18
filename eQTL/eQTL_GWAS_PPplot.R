library(magrittr)

library(dplyr)

library(ggplot2)

library(clusterProfiler)

library(org.Hs.eg.db)

library("ggrepel")

gtex=read.table("eQTL_GWAS_protein_coding.txt",sep="\t",header=T)

id_split <- gtex$variant_id |> strsplit(split="_") |> data.frame() |> t() |> as.data.frame()

gtex$hg38id <- paste(gtex[,2],id_split[,3],id_split[,4],sep = "_")

gwas_0.0001 <- read.table("assoc_results.assoc.fisher0.0001.txt",header = T, sep = "\t")

gwas_0.0001$hg38id <- paste("chr",gwas_0.0001[,1],":",gwas_0.0001[,3],"_",gwas_0.0001[,7],"_",gwas_0.0001[,4],sep="")

gwas_eQTL_merge <- merge(gwas_0.0001,gtex,by.x = "hg38id",by.y = "hg38id")

gwas_eQTL_merge <- gwas_eQTL_merge[order(gwas_eQTL_merge$pval_nominal),]

TP53BP2 <- gwas_eQTL_merge[which(gwas_eQTL_merge$Gene.name=="TP53BP2"),]



gwas_eQTL_merge <- gwas_eQTL_merge |> mutate(abs_clour=case_when(slope > 0 ~ "red",TRUE ~ "blue"))

gwas_eQTL_merge <- gwas_eQTL_merge |> mutate(highlight=case_when(abs(slope) > 0.15 ~ paste(SNP,":",Gene.name,sep=""),TRUE ~ ""))

write.table(gwas_eQTL_merge,"gwas_eQTL_merge.txt",row.names = F,quote = F)

pdf("eQTL_GWAS_PPplot.pdf",width=7,height=6)
ggplot(gwas_eQTL_merge, aes(x = -log10(P), y = -log10(pval_nominal))) +
               geom_point(data = gwas_eQTL_merge[which(gwas_eQTL_merge$Gene.name != "TP53BP2"),],
                          aes(color = I("skyblue"),
                              size = abs(slope), 
                              alpha = 1))+
               geom_point(data = gwas_eQTL_merge[which(gwas_eQTL_merge$Gene.name == "TP53BP2"),],
                          aes(x = -log10(P),
                              y = -log10(pval_nominal),
                              color=I("red"),
                              size=abs(slope),
                              alpha = 1))+
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
dev.off()