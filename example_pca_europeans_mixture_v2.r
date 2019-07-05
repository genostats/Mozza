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
# N <- 200

# x <- mix.pop(N, H, TSI = seq(0, 1, by = 0.1), IBS = seq(1, 0, by = -0.1))
# x <- select.snps(x, maf > .01)

pp <- triangle.tiling(21)
x <- mix.pop(5, H, TSI = pp$a, IBS = pp$b, CEU = pp$c)
x <- select.snps(x, maf > .01)

x.thin  <- LD.thin(x, 0.1)
K <- GRM(x.thin)
eiK <- eigen(K)
# plot(eiK$vectors, col = hsv(x@ped$TSI*0.8,1,1))
plot(eiK$vectors, col = rgb(x@ped$IBS, x@ped$TSI, x@ped$CEU))

x@ped$pheno <- rbinom( nrow(x), 1, 0.1 + 0.4*x@ped$TSI )

a <- association.test(x, method = "lm", response = "binary")
qqplot.pvalues(a)

b <- association.test(x, method = "lm", response = "binary", eigenK = eiK, p = 1)
qqplot.pvalues(b)





pp <- triangle.tiling(21)
x <- mix.pop.nuclear.families(2, 2, H, TSI = pp$a, IBS = pp$b, CEU = pp$c)
x <- select.snps(x, maf > .01)

x.thin  <- LD.thin(x, 0.1)
K <- GRM(x.thin)
eiK <- eigen(K)

expit <- function(x) 1/(1+exp(-x))
x@ped$pheno <- rbinom( nrow(x), 1, expit(  -0.25 + x@ped$TSI + rnorm.fam(x)))

a <- association.test(x, method = "lm", response = "binary")
b <- association.test(x, method = "lm", response = "binary", eigenK = eiK, p = 2)

par(mfrow = c(1,3))
plot(eiK$vectors, col = rgb(x@ped$IBS, x@ped$TSI, x@ped$CEU))
qqplot.pvalues(a, pch = ".")
qqplot.pvalues(b, pch = ".")

