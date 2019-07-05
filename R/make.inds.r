# nb.inds = nbre d'indidividus (non apparentés) à générer
# proba.haplos = un vecteur de probas pour chacun des haplotypes
# si rien n'est donné, on les prend équiprobables
# haplos = bed matrix d'haplotypes
make.inds <- function(nb.inds, haplos, proba.haplos, tile.length = 20, kinship = FALSE) {

  if(all(haplos@snps$dist == 0))
    stop("Set genetic distance between markers with set.dist !")

  if(missing(proba.haplos))
    L <- make_inds(nb.inds, tile.length, haplos@bed, haplos@snps$chr, haplos@snps$dist, kinship)
  else {
    if(is.vector(proba.haplos))
      proba.haplos <- matrix(proba.haplos, ncol = 1)
    nb.inds <- rep_len(nb.inds, ncol(proba.haplos))
    L <- make_inds_probs(nb.inds, proba.haplos, tile.length, haplos@bed, haplos@snps$chr, haplos@snps$dist, kinship)
  }

  ped <- data.frame(famid = 1:sum(nb.inds), id = 1:sum(nb.inds), father = NA, 
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

