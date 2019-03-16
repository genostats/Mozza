
# nb.inds = nbre d'indidividus (non apparentés) à générer
# haplos = bed matrix d'haplotypes
make.nuclear.families <- function(nb.fams, nb.offsprings, haplos, tile.length = 20, kinship = FALSE) {
  if(all(haplos@snps$dist == 0))
    stop("Set genetic distance between markers with set.dist !")
  
  L <- nuclear_families(nb.fams, nb.offsprings, tile.length, haplos@bed, haplos@snps$chr, haplos@snps$dist, kinship)
  
  famid <- rep( 1:nb.fams, each = (2 + nb.offsprings) )
  ids <- rep( c(1:2, 2+(1:nb.offsprings)), nb.fams )
  father <- rep(c(NA, NA, rep(1, nb.offsprings)), nb.fams)
  mother <- rep(c(NA, NA, rep(2, nb.offsprings)), nb.fams)
  sex <- rep(c(1,2, rep(NA, nb.offsprings)), nb.fams)
  ped <- data.frame(famid = famid, ids = ids, father = father, mother = mother, sex = sex, 
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

