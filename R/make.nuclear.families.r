
# nb.inds = nbre d'indidividus (non apparentés) à générer
# proba.haplos = un vecteur de probas pour chacun des haplotypes
# si rien n'est donné, on les prend équiprobables
# haplos = bed matrix d'haplotypes
make.nuclear.families <- function(nb.fams, nb.offsprings, haplos, proba.haplos, tile.length = 20, kinship = FALSE) {

  if(all(haplos@snps$dist == 0))
    stop("Set genetic distance between markers with set.dist !")

  if(missing(proba.haplos)) 
    L <- nuclear_families(nb.fams, nb.offsprings, tile.length, haplos@bed, haplos@snps$chr, haplos@snps$dist, kinship)
  else {
    if(is.vector(proba.haplos))
      proba.haplos <- matrix(proba.haplos, ncol = 1)
    nb.fams <- rep_len(nb.fams, ncol(proba.haplos))
    L <- nuclear_families_probs(nb.fams, nb.offsprings, proba.haplos, tile.length, haplos@bed, haplos@snps$chr, haplos@snps$dist, kinship)
  }

  N <- sum(nb.fams)
  famid <- rep( 1:N, each = (2 + nb.offsprings) )
  id <- rep( c(1:2, 2+(1:nb.offsprings)), N )
  father <- rep(c(NA, NA, rep(1, nb.offsprings)), N)
  mother <- rep(c(NA, NA, rep(2, nb.offsprings)), N)
  sex <- rep(c(1,2, rep(NA, nb.offsprings)), N)
  ped <- data.frame(famid = famid, id = id, father = father, mother = mother, sex = sex, 
                    pheno = NA, stringsAsFactors = FALSE)
  
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

