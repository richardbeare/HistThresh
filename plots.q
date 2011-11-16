## R script to generate plots for paper

#cthead.hist <- read.csv('hist.csv', header=F)
brainpd.hist <- read.table('brain_pd.hist', header=F)

tt <- strsplit(readLines('thresholds.txt'), ':')

ttnames <- sapply(tt, function(l){l[2]})
ttnames <- gsub(" threshold", "", ttnames)
ttnames <- gsub("[[:space:]]", "", ttnames)
ttvals <- sapply(tt, function(l){as.numeric(l[3])})
names(ttvals)=ttnames
n1 <- c("KittlerIllingworth", "Intermode", "Minimum", "IsoData", "Moments", "Triangle")
n2 <- c("Huang", "Li", "MaxEntropy", "RenyiEntropy", "Yen", "Shanbhag")

n1vals=ttvals[n1]
n2vals=ttvals[n2]

pdf("hist_resultsA.pdf")
#plot(cthead.hist, type='l', xlab="Brightness", ylab="Voxel count")
plot(brainpd.hist, type='l', xlab="Brightness", ylab="Voxel count")

for (i in 1:length(n1))
{
  abline(v=n1vals[i], col=i)
}

legend('topright', legend=n1, col=1:length(n1), lty=1)

dev.off()

pdf("hist_resultsB.pdf")
#plot(cthead.hist, type='l', xlab="Brightness", ylab="Voxel count")
plot(brainpd.hist, type='l', xlab="Brightness", ylab="Voxel count")

for (i in 1:length(n2))
{
  abline(v=n2vals[i], col=i)
}

legend('topright', legend=n2, col=1:length(n2), lty=1)

dev.off()


