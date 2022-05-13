# Author: Preshita 
# Description: For each run, plot the cumulative distribution and retrieve whitelist

setwd('/projectnb/bf528/students/preshita/project-4/data/')
counts <- read.table("merged_counts.txt", row.names = 'V2')

colnames(counts) <- 'count'

# Order by counts for the cdf plot
counts <- counts[order(counts$count),]

#calculate empirical CDF of data
plot(ecdf(counts$count), main='CDF of Barcode Frequency',xlim=c(0,2000))
                                                              

plot(ecdf(counts$count), main='CDF of Barcode Frequency',
     xlab='Rank')

writeLines(as.character(rownames(counts >= 100)), file('whitelist.txt'), sep="\n")
