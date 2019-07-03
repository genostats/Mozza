
mix.pop <- function(nb.inds, haplos, population = haplos@ped$population, tile.length = 20, kinship = FALSE, ...) {
  probs <- list(...)
  p <- table(population)
  proba.haplos <- numeric( nrow(haplos) )
  for(x in names(p)) {
    if(is.null(probs[[x]]))
      proba.haplos[ population == x ] <- 0
    else
       proba.haplos[ population == x ] <- probs[[x]]/p[[x]]
  }
  make.inds(nb.inds, haplos, proba.haplos, tile.length, kinship)
}
