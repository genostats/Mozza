require(Mozza)

## Charger les haplotypes 'templates'
H <- read.bed.matrix("~/PROJECTS/PCA/2019/1KG_haplos")

## Exemple 1, segment : mélange continu de deux populations 
## Ici 11 sous-populations dans les proportions données par les vecteurs TSI / IBS
## Avec chaque fois, 100 individus 
N <- 100
x <- mix.pop(N, H, TSI = seq(0, 1, by = 0.1), IBS = seq(1, 0, by = -0.1))
x <- select.snps(x, maf > 0.05)

# une ACP du rouge (TSI) au vert (IBS)
# la matrice x a dans son @ped les proportions de TSI / IBS de chaque individu
# on les utilise pour la couleur.
x.thin <- LD.thin(x, .1)
K <- GRM(x.thin)
eigenK <- eigen(K)
plot(eigenK$vectors, col = rgb(x@ped$TSI, x@ped$IBS, 0), pch = 16)


## Exemple 2, triangle : mélange de trois populations
## Les sous-populations sont distribuées sur un triangle équilatéral
## avec 21 sous pops sur les cotés. Le premier plot permet de visualiser
## la répartition des 231 = 21 * 22 / 2 populations résultantes. La
## proportion du mélange dépend des coordonnées barycentriques des 
## sous pops.
pp <- triangle.tiling(21)
plot( triangle.plot(pp$a, pp$b, pp$c), asp = 1 )

# pour aller vite, 5 individus dans chaque sous population -> 5*231 = 1155 individus
x <- mix.pop(5, H, TSI = pp$a, IBS = pp$b, CEU = pp$c)
x <- select.snps(x, maf > 0.05)

# ACP, meme démarche
x.thin <- LD.thin(x, .1)
K <- GRM(x.thin)
eigenK <- eigen(K)
plot(eigenK$vectors, col = rgb(x@ped$TSI, x@ped$IBS, x@ped$CEU), pch = 16)


## Exemple 3, carré : 4 populations (il serait très intéressant de faire un 
## tétraèdre plutôt qu'un carré mais j'ai essayé de faire quelque chose de 
## plan... choix discutable)
f <- function(x,y) {
  c( (1-y)*c(x, 1-x), y*c(x, 1-x) )
}

t <- rep(seq(0,1,length=11), 11)
u <- rep(seq(0,1,length=11), each = 11)
pp <- t(mapply(f, t, u))
x <- mix.pop(10, H, TSI = pp[,1], IBS = pp[,2], CEU = pp[,3], FIN = pp[,4])
x <- select.snps(x, maf > 0.05)

x.thin <- LD.thin(x, .1)
K <- GRM(x.thin)
eigenK <- eigen(K)
plot(eigenK$vectors, col = rgb(x@ped$TSI, x@ped$IBS, x@ped$CEU), pch = 16)

