# nb.inds = nbre d'indidividus (non apparentés) à générer
# haplos = bed matrix d'haplotypes
make.inds <- function(nb.inds, haplos, tile.length = 20, kinship = FALSE) {
  if(all(haplos@snps$dist == 0))
    stop("Set genetic distance between markers with set.dist !")
  
  L <- make_inds(nb.inds, tile.length, haplos@bed, haplos@snps$chr, haplos@snps$dist, kinship)
  
  ped <- data.frame(famid = 1:nb.inds, ids = 1:nb.inds, father = NA, 
                    mother = NA, sex = NA, pheno = NA, stringsAsFactors = FALSE)
  
  x <- new("bed.matrix", bed = L$bed, snps = haplos@snps, ped = ped, p = NULL, mu = NULL,
           sigma = NULL, standardize_p = FALSE, standardize_mu_sigma = FALSE )
  if(getOption("gaston.auto.set.stats", TRUE))
    x <- set.stats(x)
  if(kinship) {
    L$bed = x
    return(L);
  }
  x
}

