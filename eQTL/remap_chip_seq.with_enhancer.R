library(magrittr)

ReMap <- read.table("remap_chip_seq.txt",header=T,sep="\t")

ReMap$enhancer <- "enhancer:"
enhancer <- read.table("enhancer_position.txt",header=T,sep="\t")


for(i in c(1:nrow(ReMap))){
  ReMap_line <- c(ReMap$chromStart[i]:ReMap$chromEnd[i])
  for (j in c(1:nrow(enhancer))) {
    enhancer_line <- c(enhancer$start[j]:enhancer$end[j])
    inter <- intersect(ReMap_line,enhancer_line) %>% length()
    if(inter>0){ReMap$enhancer[i] <- paste(ReMap$enhancer[i],enhancer$name[j],sep=",")}
  }
}


write.table(ReMap,"remap_chip_seq.with_enhancer.txt",col.names = T, row.names = F, quote = F, sep = "\t")