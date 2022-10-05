library(edgeR)
library(org.Hs.eg.db)
library(VennDiagram)
library(dbplyr)
library(stringr)
library(pheatmap)
library(gridExtra)
library(factoextra)
library(tidyverse)

foldChange=log2(1.5)
padj=0.05
###data prepaired

#HepG2

rt <- read.table("symbol.protein_coding.hepg2.txt",
                    sep="\t",header=T,check.names=F)

species="hsa"
org.db="org.Hs.eg.db"
group1=c(rep("groupA",3))
group2=c(rep("groupB",3))
group=c(group1,group2)
design <- model.matrix(~group)

rt <-aggregate(rt[,-1],by=list(rt[,1]),FUN=mean)
rownames(rt)=rt[,1]
rt=rt[,-1]

all <- rt

#HepG2 siNC_siTP53BP2

rt <- all[,c(1,2,3,7,8,9)]

rt=as.matrix(rt)
exp=rt
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)

data=avereps(data)
data=data[rowMeans(data)>1,]

y <- DGEList(counts=data,group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)

et <- exactTest(y,pair = c("groupA","groupB"))

ordered_tags <- topTags(et, n=100000)
allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
diff=allDiff

diff_g2_siNC_siTP53BP2 <- diff[,c(1,3,4)]

diff_g2_siNC_siTP53BP2[,4] <- rownames(diff_g2_siNC_siTP53BP2)

diffSig_g2_siNC_siTP53BP2 <-  diff[(diff$FDR < padj & (diff$logFC>foldChange | 
                                                       diff$logFC<(-foldChange))),]

write.csv(diffSig_g2_siNC_siTP53BP2,"diffSig_g2_siNC_siTP53BP2.csv")

#HepG2 siNC_siTP53BP2_1h

rt <- all[setdiff(rownames(all),rownames(diffSig_g2_siNC_siTP53BP2)),c(4,5,6,10,11,12)]

rt=as.matrix(rt)
exp=rt
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)

data=avereps(data)
data=data[rowMeans(data)>1,]

y <- DGEList(counts=data,group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)

et <- exactTest(y,pair = c("groupA","groupB"))

ordered_tags <- topTags(et, n=100000)
allDiff <- ordered_tags$table
allDiff <- allDiff[is.na(allDiff$FDR)==FALSE,]
diff <- allDiff

diff_g2_siNC_siTP53BP2_1h <- diff[,c(1,3,4)]

diff_g2_siNC_siTP53BP2_1h[,4] <- rownames(diff_g2_siNC_siTP53BP2_1h)



diffsig_g2_siNC_siTP53BP2_1h <- diff[(diff$FDR < padj & (diff$logFC>foldChange | 
                                                           diff$logFC<(-foldChange))),]

write.csv(diffsig_g2_siNC_siTP53BP2_1h,"diffsig_g2_siNC_siTP53BP2_1h.csv")
													   
##Heatmap for Type I interferon pathway genes

newData=y$pseudo.counts

Type_I_interferon_pathway_genes <- c(read.table("Type_I_interferon_pathway_genes.txt",header = F, sep = "\t")[,1])

df <- newData[Type_I_interferon_pathway_genes,]

col <- colorRampPalette(c("blue", "white", "red"))(256)
patient_class=RowSideColors =  c(rep("purple", 3), rep("orange", 3))

pdf("Type_I_interferon_pathway_genes.pdf")
heatmap(df[rev(rownames(df)),],
        col = col,
        Rowv = NA,
        Colv = NA,
        ColSideColors = patient_class)
dev.off()													   
													   
diffsig_Type_I_interferon_pathway_genes <- diffsig_g2_siNC_siTP53BP2_1h[Type_I_interferon_pathway_genes,]								   
write.csv(diffsig_Type_I_interferon_pathway_genes,"diffsig_Type_I_interferon_pathway_genes.csv")

#HepAD38
rt <- read.table("symbol.protein_coding.hepad38.txt",
                 sep="\t",header=T,check.names=F)

rt <-aggregate(rt[,-1],by=list(rt[,1]),FUN=mean)
rownames(rt)=rt[,1]
rt=rt[,-1]

all_ad38 <- rt

#HepAD38 siNC_siTP53BP2
rt <- all_ad38[,c(1,5,6,7,8,9)]

rt=as.matrix(rt)
exp=rt
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)

data=avereps(data)
data=data[rowMeans(data)>1,]

y <- DGEList(counts=data,group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)

et <- exactTest(y,pair = c("groupA","groupB"))

ordered_tags <- topTags(et, n=100000)
allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
diff=allDiff

diff_ad38_siNC_siTP53BP2 <- diff[,c(1,3,4)]

diff_ad38_siNC_siTP53BP2[,4] <- rownames(diff_ad38_siNC_siTP53BP2)

diffSig_ad38_siNC_siTP53BP2 = diff[(diff$FDR < padj & (diff$logFC>foldChange | 
                                                          diff$logFC<(-foldChange))),]
#Volcano plot (Sup fig3A)
pdf(file="vol.pdf")
allDiff$FDR[-log10(allDiff$FDR) >= 100] <- 10^-100
xMax=max(-log10(allDiff$FDR))+1
yMax=12
plot(-log10(allDiff$FDR), allDiff$logFC, xlab="-log10(FDR)",ylab="logFC",
     main="Volcano", xlim=c(0,xMax),ylim=c(-yMax,yMax),yaxs="i",pch=20, cex=0.4)
diffSub=allDiff[allDiff$FDR<padj & allDiff$logFC>foldChange,]
points(-log10(diffSub$FDR), diffSub$logFC, pch=20, col="red",cex=0.4)
diffSub=allDiff[allDiff$FDR<padj & allDiff$logFC<(-foldChange),]
points(-log10(diffSub$FDR), diffSub$logFC, pch=20, col="green",cex=0.4)
abline(h=0,lty=2,lwd=3)
dev.off()
											  														  

#LM3 shNC_shTP53BP2

rt <- read.table("lm3_edgerOut.tsv",
                 sep="\t",header=T,check.names=F)

diff_lm3_siNC_siTP53BP2 <- rt[,c(1,2,5,6)]

rownames(diff_lm3_siNC_siTP53BP2) <- diff_lm3_siNC_siTP53BP2[,1]

diff_lm3_siNC_siTP53BP2 <- diff_lm3_siNC_siTP53BP2[,-1]

colnames(diff_lm3_siNC_siTP53BP2)[2] <- "PValue"

diff_lm3_siNC_siTP53BP2[,4] <- rownames(diff_lm3_siNC_siTP53BP2)

rt <- rt[which(rt$logFC>=foldChange | rt$logFC<(-foldChange)),]

rt <- rt[which(rt$FDR<=padj),]


diffSig_lm3_siNC_siTP53BP2 <- rt[,-1]
rownames(diffSig_lm3_siNC_siTP53BP2) <- rt[,1]

all_lm3 <- read.table("lm3.normalized.csv",
                      sep=",",header=T,check.names=F)
rownames(all_lm3) <- all_lm3[,1] 
all_lm3 <- all_lm3[,-1]


###venn plot

gene_logFC <- function(diff,foldChange,TP53BP2_DOWN){
  diff$change <- "UPUP"
  diff$change[which(diff$logFC< (-foldChange))] <- "DOWN"
  name <- paste(row.names(diff),diff$change,sep="_")
  return(name)
}



HepG2_1h <- gene_logFC(diffSig_g2_siNC_siTP53BP2_1h,foldChange) %>% setdiff(.,y="TP53BP2_DOWN")

HepG2_0h <- gene_logFC(diffSig_g2_siNC_siTP53BP2,foldChange) %>% setdiff(.,y="TP53BP2_DOWN")

HepAD38 <- gene_logFC(diffSig_ad38_siNC_siTP53BP2,foldChange) %>% setdiff(.,y="TP53BP2_DOWN")

LM3 <- gene_logFC(diffSig_lm3_siNC_siTP53BP2,foldChange) %>% setdiff(.,y="TP53BP2_DOWN")




mydata <- list(
  HepG2_1h=HepG2_1h, 
  HepG2_0h=HepG2_0h,  
  HepAD38=HepAD38,
  LM3=LM3 
)


#venn plot (Figure 7A)

venn.plot <- venn.diagram(
  mydata,
  filename = NULL,
  fill = c("#0073C2FF","#E69F00","#009E73","#CD534CFF"),
  alpha = 0.30,
  lwd = 1,
  cat.col = c("steelblue", "orange", "darkgreen", "tomato"),
  cat.cex = 1,
  cat.fontface = "bold",
  cat.pos = 0,
  cat.dist = 0.03,
  cex = 0.8,
  fontface = "italic",
  margin = 0.2
)


pdf(file="venn.pdf")
grid.draw(venn.plot)
dev.off()

###heatmap for intersect different expression genes


merge_gene <- intersect(
    intersect(HepG2_1h,HepG2_0h),
    intersect(HepAD38,LM3))


merge_gene2 <- c("SOCS2","GLIPR1","TJP3")

#Build expression lists

g2 <- read.table("normalizeExp.g2.txt",
                       sep="\t",header=T,check.names=F)

rownames(g2) <- g2[,1]
g2 <- g2[,-1]


##PCA plot (Figure 4A)
g2 <- g2[apply(g2, 1, var)!=0,]
mads <- apply(g2, 1, mad)
g2 <- g2[rev(order(mads)),]
g2_t <- t(g2)
variableL <- ncol(g2_t)
pca <- prcomp(g2_t[,1:variableL], scale=T)
names<-colnames(g2)
pdf("pca.pdf",width = 1000,height = 1000)
fviz_pca_ind(pca, col.ind=names, mean.point=F, addEllipses = T, legend.title="Groups")
dev.off()





ad38 <- read.table("normalizeExp.ad38.txt",
                   sep="\t",header=T,check.names=F)

rownames(ad38) <- ad38[,1]
ad38 <- ad38[,-1]

g2 <- g2[merge_gene2,]

ad38_siNC_siTP53BP2 <- ad38[merge_gene2,]

lm3_siNC_siTP53BP2 <- all_lm3[merge_gene2,]

#Heatmap (Figure 7B)

#g2_heatmap
df <- as.matrix(g2) 
col <- colorRampPalette(c("#16499D", "#FBFCFE", "#E72119"))(256)
annotation_col<-data.frame(siTP53BP2=factor(c(rep("siControl",6),rep("siTP53BP2",6))),
                           IFN=factor(c(rep("Control",3),rep("INF-¦Á",3),
                                        rep("Control",3),rep("INF-¦Á",3)))
                           )
rownames(annotation_col)<-colnames(df)


g2_heatmap <- pheatmap(df,
                      legend = F,
                      scale = "row",
                      color = col,
                      cluster_rows = F,
                      cluster_cols = F,
                      annotation_col = annotation_col,
                      show_colnames = F)

#ad38_heatmap
df <- as.matrix(ad38_siNC_siTP53BP2) 
col <- colorRampPalette(c("#16499D", "#FBFCFE", "#E72119"))(256)
annotation_col<-data.frame(siTP53BP2=factor(c(rep("siControl",3),rep("siTP53BP2",3)))
)
rownames(annotation_col)<-colnames(df)


ad38_heatmap <- pheatmap(df,
                       legend = F,
                       scale = "row",
                       color = col,
                       cluster_rows = F,
                       cluster_cols = F,
                       annotation_col = annotation_col,
                       show_colnames = F)


#lm3_heat_map

df <- as.matrix(lm3_siNC_siTP53BP2) 
col <- colorRampPalette(c("#16499D", "#FBFCFE", "#E72119"))(256)
annotation_col<-data.frame(siTP53BP2=factor(c(rep("siControl",3),rep("siTP53BP2",3)))
)
rownames(annotation_col)<-colnames(df)


lm3_heatmap <- pheatmap(df,
                         legend = F,
                         scale = "row",
                         color = col,
                         cluster_rows = F,
                         cluster_cols = F,
                         annotation_col = annotation_col,
                         show_colnames = F)

dev.off()

#output heatmap 
pdf("heatmap.pdf")
plot_list=list()
grid.arrange(arrangeGrob(g2_heatmap[[4]]), 
             arrangeGrob(ad38_heatmap[[4]],lm3_heatmap[[4]], ncol=2),
             ncol=1)
dev.off()


#merge Volcano

changeName <- function(data,name){
  rownames(data) <- paste0(rownames(data),"_",name)
  return(data)
  }


diff_ad38_siNC_siTP53BP2[,5] <- "HepAD38"
diff_g2_siNC_siTP53BP2[,5] <- "HepG2"
diff_g2_siNC_siTP53BP2_1h[,5] <- "HepG2-IFN"

diff_ad38_siNC_siTP53BP2 <- changeName(diff_ad38_siNC_siTP53BP2,"HepAD38")
diff_g2_siNC_siTP53BP2 <- changeName(diff_g2_siNC_siTP53BP2,"HepG2")
diff_g2_siNC_siTP53BP2_1h <- changeName(diff_g2_siNC_siTP53BP2_1h,"HepG2-IFN")



merge.diff <- rbind(diff_ad38_siNC_siTP53BP2,
                    diff_g2_siNC_siTP53BP2,
                    diff_g2_siNC_siTP53BP2_1h)

merge.diff[,6] <- "Stable"
merge.diff[which(merge.diff$logFC > foldChange & merge.diff$FDR < padj),6] <- "Up"
merge.diff[which(merge.diff$logFC < (-foldChange) & merge.diff$FDR < padj),6] <- "Down"


colnames(merge.diff)[c(4,5,6)] <- c("geneids","contrast","change")

my_genes <- Type_I_interferon_pathway_genes

ggplot(filter(merge.diff,!geneids %in% my_genes)) +
  geom_point(aes(contrast, logFC,color = change, size = -log10(FDR), alpha = -log10(FDR)),
             position = "jitter")+
  geom_point(data = filter(merge.diff,geneids %in% my_genes),
             aes(contrast, logFC), 
             position = "jitter",
             size = 4,
             shape = 21,
             color = "black",
             fill = "#fee08b")+
  geom_hline(yintercept = c(-foldChange, foldChange), 
             linetype = 'dotdash',
             color = 'grey30')+
  scale_color_manual(values = c(Up = "#e6550d", Down = "#3182bd",Stable="#636363")) +
  scale_alpha(range = c(0.3, 1)) +
  labs(x = 'Group', y = 'Log2(fold change)') +
  theme_bw()
