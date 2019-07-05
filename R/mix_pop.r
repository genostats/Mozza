
make.proba.haplos <- function(population, probs) {
  N <- max(sapply(probs, length))
  if( N == 0 )
    stop("Please give proportions vector for some of the populations")
  probs <- lapply( probs, rep_len, N)

  p <- table(population)
  proba.haplos <- matrix(0.0, nrow = length(population), ncol = N)
  for(x in names(probs)) {
    for(j in 1:N) 
      proba.haplos[ population == x , j] <- probs[[x]][j]/p[[x]]
  }
  proba.haplos
}

mix.pop <- function(nb.inds, haplos, population = haplos@ped$population, tile.length = 20, kinship = FALSE, ...) {
  if(length(population) != nrow(haplos))
    stop("Dimensions mismatch")
  probs <- list(...)
  proba.haplos <- make.proba.haplos(population, probs)
  nb.inds <- rep_len(nb.inds, ncol(proba.haplos))
  x <- make.inds(nb.inds, haplos, proba.haplos, tile.length, kinship)
  for(a in names(probs)) {
    x@ped[[a]] <- as.vector(mapply(rep, probs[[a]], nb.inds))
  }
  x
}

mix.pop.nuclear.families <- function(nb.fams, nb.offsprings, haplos, population = haplos@ped$population, tile.length = 20, kinship = FALSE, ...) {
  if(length(population) != nrow(haplos))
    stop("Dimensions mismatch")
  probs <- list(...)
  proba.haplos <- make.proba.haplos(population, probs)
  nb.fams <- rep_len(nb.fams, ncol(proba.haplos))
  x <- make.nuclear.families(nb.fams, nb.offsprings, haplos, proba.haplos, tile.length, kinship)
  for(a in names(probs)) {
    x@ped[[a]] <- as.vector(mapply(rep, probs[[a]], nb.fams*(2+nb.offsprings)))
  }
  x
}

