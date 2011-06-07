## R script to generate plots for paper

cthead.hist <- read.csv('hist.csv', header=F)

tt <- strsplit(readLines('thresholds.txt'), ':')

ttnames <- sapply(tt, function(l){l[2]})
ttnames <- gsub(" threshold", "", ttnames)
ttvals <- sapply(tt, function(l){as.numeric(l[3])})

pdf("hist_results.pdf")
plot(cthead.hist, type='l', xlab="Brightness", ylab="Voxel count")

for (i in 1:length(ttnames))
{
  abline(v=ttvals[i], col=i)
}

legend('topright', legend=ttnames, col=1:length(ttnames), lty=1)
dev.off()

