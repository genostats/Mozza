mix.pop.nuclear.families <- function(nb.fams, nb.offsprings, haplos, population = haplos@ped$population, tile.length = 20, 
                                     kinship = FALSE, fraternity = FALSE, ...) {
  if(length(population) != nrow(haplos))
    stop("Dimensions mismatch")
  probs <- list(...)
  proba.haplos <- make.proba.haplos(population, probs)
  nb.fams <- rep_len(nb.fams, ncol(proba.haplos))
  x <- make.nuclear.families(nb.fams, nb.offsprings, haplos, proba.haplos, tile.length, kinship, fraternity)
  for(a in names(probs)) {
    x$bed@ped[[a]] <- as.vector(mapply(rep, probs[[a]], nb.fams*(2+nb.offsprings)))
  }
  x
}

