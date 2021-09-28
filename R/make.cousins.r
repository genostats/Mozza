
# nb.inds = nbre d'indidividus (non apparentés) à générer
# haplos = bed matrix d'haplotypes
make.cousins <- function(n, haplos, tile.length = 20, kinship = FALSE, fraternity = FALSE) {
  if(all(haplos@snps$dist == 0))
    stop("Set genetic distance between markers with set.dist !")
  
  L <- cousins_1stdegree(n, tile.length, haplos@bed, haplos@snps$chr, haplos@snps$dist, kinship, fraternity)
  
  famid <- rep( 1:n, each = 2)
  id <- rep( 1:2, n )
  ped <- data.frame(famid = famid, id = id, father = NA, mother = NA, sex = NA, 
                    pheno = NA, stringsAsFactors = FALSE)
  
  x <- new("bed.matrix", bed = L$bed, snps = haplos@snps, ped = ped, p = NULL, mu = NULL,
           sigma = NULL, standardize_p = FALSE, standardize_mu_sigma = FALSE )
  
  if(getOption("gaston.auto.set.stats", TRUE))
    x <- set.stats(x)
  
  L$bed <- x
  L
}

