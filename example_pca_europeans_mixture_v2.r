require(Mozza)
require(gaston.utils)

H <- read.bed.matrix("~/PROJECTS/PCA/2019/1KG_haplos")
H <- set.dist(H, HumanGeneticMap::genetic.map.b37)

# si les pops n'ont pas été bien mises dans 1KG_haplos.ped... :/
if( !all( H@ped$famid == substr(H@ped$famid,1,7) ) ) {
  H@ped$famid <- substr(H@ped$famid,1,7)
  m <- match(H@ped$famid, gaston.utils::KG.samples$sample)
  H@ped$population       <- droplevels( gaston.utils::KG.samples$population[m] )
  H@ped$super.population <- droplevels( gaston.utils::KG.samples$super.population[m] )
}

# N <- 400
N <- 200

x <- mix.pop(N, H, TSI = 0, IBS = 1)
x@ped$TSI <- 0

for(p in seq(0.1,1,by = 0.1)) {
  # x1 <-  mix.pop(N,  H, TSI = p, IBS = 1-p)
  x1 <-  mix.pop(N,  H, TSI = p, IBS = 1-p, CEU = p*(1-p))
  x1@ped$TSI <- p
  x <- rbind(x, x1)
}

x@ped$pheno <- rbinom( nrow(x), 1, 0.1 + 0.4*x@ped$TSI )

a <- association.test(x, method = "lm", response = "binary")
qqplot.pvalues(a)

x.thin  <- LD.thin(x, 0.1)
K <- GRM(x.thin)
eiK <- eigen(K)

plot(eiK$vectors, col = hsv(x@ped$TSI*0.8,1,1))

b <- association.test(x, method = "lm", response = "binary", eigenK = eiK, p = 1)
qqplot.pvalues(b)

