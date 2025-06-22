data <- read.table(file="own.mds", header=TRUE)
 race <- read.table(file="raceown.txt", header=FALSE)
colnames(race)=c("FID","IID","race")
datafile <- merge(data, race, by=c("IID"))
datafile <- datafile[,-14]
head(datafile)

# PC1 vs PC2 plot
pdf("MDS.own.PC1.PC2.pdf", width=7, height=7)
plot(NA, xlim=c(-0.05,0.05), ylim=c(-0.05,0.05), xlab="PC 1", ylab="PC 2")
points(datafile[datafile$race == 1, 4], datafile[datafile$race == 1, 5],
       pch=1, cex=0.5, col="blue")
points(datafile[datafile$race == 2, 4], datafile[datafile$race == 2, 5],
       pch=3, cex=0.5, col="red")
abline(v=-0.035, lty=3)
abline(h=0.035, lty=3)
legend("topright", pch=c(1,3), c("control","case"), col=c("blue","red"), bty="o", cex=1)
dev.off()

# PC2 vs PC3 plot
pdf("MDS.own.PC2.PC3.pdf", width=7, height=7)
plot(NA, xlim=c(-0.05,0.05), ylim=c(-0.05,0.05), xlab="PC 2", ylab="PC 3")
points(datafile[datafile$race == 1, 5], datafile[datafile$race == 1, 6],
       pch=1, cex=0.5, col="blue")
points(datafile[datafile$race == 2, 5], datafile[datafile$race == 2, 6],
       pch=3, cex=0.5, col="red")
abline(v=-0.035, lty=3)
abline(h=0.035, lty=3)
legend("topright", pch=c(1,3), c("control","case"), col=c("blue","red"), bty="o", cex=1)
dev.off()
