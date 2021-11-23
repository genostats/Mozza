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
