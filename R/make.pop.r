
# n0 = nb individus en génération 0
# n.gen = nb générations à simuler
# n.keep = nb de générations à garder
# lambda = paramètre loi de Poisson (nb enfants par couple) (lambda = 2 suggéré)

# proba.haplos = un vecteur de probas pour chacun des haplotypes
# si rien n'est donné, on les prend équiprobables
# haplos = bed matrix d'haplotypes
make.pop <- function(n0, n.gen, n.keep, lambda, haplos, proba.haplos, tile.length = 20, 
                                  kinship = FALSE, fraternity = FALSE) {

  if(all(haplos@snps$dist == 0))
    stop("Set genetic distance between markers with set.dist !")
  if(n.keep > n.gen)
    stop("n.gen must be greater or equal to n.keep")

  if(missing(proba.haplos)) 
    L <- population(n0, n.gen, n.keep, lambda, tile.length, haplos@bed, haplos@snps$chr, haplos@snps$dist, kinship, fraternity)
  else {
    if(is.vector(proba.haplos))
      proba.haplos <- matrix(proba.haplos, ncol = 1)
    stop("NOT IMPLEMENTED YET")
  }

  N <- L$N
  famid <- rep(NA_integer_,N)
  # renumbering ids from 1 to N
  id <- 1:N 
  father <- match(L$father, L$id)
  mother <- match(L$mother, L$id)
  
  sex <- rep(NA_integer_, N)
  sex[ father ] <- 1
  sex[ mother ] <- 2
  ped <- data.frame(famid = famid, id = id, father = father, mother = mother, sex = sex, 
                    pheno = NA, stringsAsFactors = FALSE)
  
  x <- new("bed.matrix", bed = L$bed, snps = haplos@snps, ped = ped, p = NULL, mu = NULL,
           sigma = NULL, standardize_p = FALSE, standardize_mu_sigma = FALSE )
  
  if(getOption("gaston.auto.set.stats", TRUE))
    x <- set.stats(x)

  L$bed <- x
  L$id <- L$father <- L$mother <- NULL
  L
}

