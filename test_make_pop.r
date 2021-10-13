require(Mozza)
## Charger les haplotypes 'templates'
H <- read.bed.matrix("~/PROJECTS/PCA/2019/1KG_haplos")

# n0 = 200, 10 générations...
L <- make.pop(200, 10, 2, 2, H, kinship = TRUE, fraternity = TRUE)
x1 <- select.snps(L$bed, maf > 0.05)
x1
K <- GRM(x1)
plot(K, L$kinship, pch = "."); abline(0,1)
D <- DM(x1)
plot(D, L$fraternity, pch = "."); abline(0,1)

