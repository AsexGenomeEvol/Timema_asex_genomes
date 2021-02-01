#!/usr/bin/Rscript
require(bio3d)
align <- read.fasta("aligntmp")                                                                                                                                                                      
dist <- seqidentity(align, normalize=TRUE) 
n <- nrow(dist)
distances <- matrix(NA, ncol=2, nrow=(n-1))
for (i in 1:(n-1)) {
	distances[i,1] <- rownames(dist)[i]
	distances[i,2] <- dist[i,n]
}
write.table(distances, file = "distances", quote = F, sep = " ", row.names = F)
