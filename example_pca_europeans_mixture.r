H <- read.bed.matrix("~/PROJECTS/PCA/2019/1KG_haplos")
H <- set.dist(H, HumanGeneticMap::genetic.map.b37)

# si les pops n'ont pas été bien mises dans 1KG_haplos.ped... :/
if( !all( H@ped$famid == substr(H@ped$famid,1,7) ) ) {
  H@ped$famid <- substr(H@ped$famid,1,7)
  m <- match(H@ped$famid, gaston.utils::KG.samples$sample)
  H@ped$population       <- droplevels( gaston.utils::KG.samples$population[m] )
  H@ped$super.population <- droplevels( gaston.utils::KG.samples$super.population[m] )
}


p0 <- ifelse( H@ped$population %in% c("TSI", "IBS"), 1, 0)
x0 <- make.inds(200, H, p0)
x0@ped$population <- "MIX"
x0@ped$famid <- paste0("MIX.", x0@ped$famid)

p1 <- ifelse( H@ped$population %in% c("TSI"), 1, 0)
x1 <- make.inds(200, H, p1)
x1@ped$population <- "TSI"
x1@ped$famid <- paste0("TSI.", x1@ped$famid)

p2 <- ifelse( H@ped$population %in% c("IBS"), 1, 0)
x2 <- make.inds(200, H, p2)
x2@ped$population <- "IBS"
x2@ped$famid <- paste0("IBS.", x2@ped$famid)

x <- rbind(x0,x1,x2)


K <- GRM(x)
eiK <- eigen(K)
plot(eiK$vectors, col = x@ped$population)
L <- bed.loadings(x, eiK$vectors[,1:2])
plot(L[,1], col = x@snps$chr, pch = ".") 
