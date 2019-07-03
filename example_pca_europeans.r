H <- read.bed.matrix("~/PROJECTS/PCA/2019/1KG_haplos")
H <- set.dist(H, HumanGeneticMap::genetic.map.b37)

#x <- make.inds(2000, H, tile.length = 40)
x <- make.inds(600, H, tile.length = 40)
K <- GRM(x)
eiK <- eigen(K)
plot(eiK$vectors)
L <- bed.loadings(x, eiK$vectors[,1:2])
plot(L[,1], col = x@snps$chr, pch = ".") # HLA domine !
