
make.pairs <- function(nb.pairs, haplos, shared.length, unshared.length, tile.length = 20) {
  if(all(haplos@snps$dist == 0))
    stop("Set genetic distance between markers with set.dist !")
  
  L <- make_pairs(nb.pairs, shared.length, unshared.length, tile.length, 
                  haplos@bed, haplos@snps$chr, haplos@snps$dist)

  ped <- data.frame(famid = rep(1:nb.pairs, each=2), id = c("A", "B"), father = NA, 
                    mother = NA, sex = NA, pheno = NA, stringsAsFactors = FALSE)
  
  x <- new("bed.matrix", bed = L$bed, snps = haplos@snps, ped = ped, p = NULL, mu = NULL,
           sigma = NULL, standardize_p = FALSE, standardize_mu_sigma = FALSE )
  
  if(getOption("gaston.auto.set.stats", TRUE))
    x <- set.stats(x)
  
  list("bed.matrix" = x, "pairs.kinship" = L$kin)
}
