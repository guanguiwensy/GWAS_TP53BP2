data<- read.table(file="MDS_merge2.mds",header=TRUE)
race<- read.table(file="racefile.txt",header=TRUE)
datafile<- merge(data,race,by=c("IID"))
datafile=datafile[,-14]
head(datafile)
pdf("MDS.own.PC1.PC2.pdf",width=7,height=7)
for (i in 1:nrow(datafile))
{
if (datafile[i,14]==1) {plot(datafile[i,4],datafile[i,5],type="p",xlim=c(-0.1,0.1),ylim=c(-0.1,0.1),xlab="PC 1",ylab="PC 2",pch=1,cex=0.5,col="blue")}
par(new=T)
if (datafile[i,14]==2) {plot(datafile[i,4],datafile[i,5],type="p",xlim=c(-0.1,0.1),ylim=c(-0.1,0.1),xlab="PC 1",ylab="PC 2",pch=3,cex=0.5,col="red")}
par(new=T)
#text(datafile[i,4],datafile[i,5],datafile[i,1])
}

abline(v=-0.035,lty=3)
abline(h=0.035,lty=3)
legend("topright", pch=c(1,3),c("control","case"),col=c("blue","red"),bty="o",cex=1)

dev.off()

pdf("MDS.own.PC2.PC3.pdf",width=7,height=7)
for (i in 1:nrow(datafile))
{
if (datafile[i,14]==1) {plot(datafile[i,5],datafile[i,6],type="p",xlim=c(-0.1,0.1),ylim=c(-0.1,0.1),xlab="PC 2",ylab="PC 3",pch=1,cex=0.5,col="blue")}
par(new=T)
if (datafile[i,14]==2) {plot(datafile[i,5],datafile[i,6],type="p",xlim=c(-0.1,0.1),ylim=c(-0.1,0.1),xlab="PC 2",ylab="PC 3",pch=3,cex=0.5,col="red")}
par(new=T)
#text(datafile[i,4],datafile[i,5],datafile[i,1])
}

abline(v=-0.035,lty=3)
abline(h=0.035,lty=3)
legend("topright", pch=c(1,3),c("control","case"),col=c("blue","red"),bty="o",cex=1)

dev.off()
