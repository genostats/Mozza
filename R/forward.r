#' Forward time 
#'
#' @param x a vector of zygotes
#' @param n.gen number of generations to simulate
#' @param n.keep number of generations to keep 
#' @param lambda mean number of offspring by pair (default is 2)
#'
#' @details Starting from a population in x, n.gen generations are simulated by random mating, each mate pair having 
#' a random number of offspring, drawn in a Poisson distribution of parameter lambda. The last n.keep generations 
#' are sent to the user.
#'
#' @return A list with named elements 'zygotes', 'father', 'mother', 'sex'.
#'
#' @examples 
#' # installs KGH if not already installed
#' if(!require("KGH")) install.packages("KGH", repos="https://genostats.github.io/R/")
#' # loads a bed matrix of 1006 european haplotypes
#' filepath <- system.file("extdata", "1KG_haplos.bed", package = "KGH")
#' H <- read.bed.matrix(filepath)
#' # generates 100 individuals
#' set.seed(1)
#' x <- make.inds(100, nrow(H))
#' # forward 10 generations, keep only last two generations
#' fo <- forward(x, 10, 2)
#' # Lots of related individuals !
#' K <- kinship.matrix(fo$zygotes)
#' round(K[1:10, 1:10], 1)
#' # to get a bed matrix, use drop.genotypes (here on SNPs 1 to 100)
#' ybed <- drop.genotypes(fo$zygotes, H[,1:100])
#' # integrate father / mother / sex information sent by the function
#' ybed@ped$father <- fo$father
#' ybed@ped$mother <- fo$mother
#' ybed@ped$sex <- fo$sex

#' @export
forward <- function(zygotes, n.gen, n.keep, lambda = 2) {

  if(n.keep > n.gen)
    stop("n.gen must be greater or equal to n.keep")

  L <- forward_(zygotes, n.gen, n.keep, lambda)

  N <- L$N
  # renumbering ids from 1 to N
  id <- 1:N 
  father <- match(L$father, L$id)
  mother <- match(L$mother, L$id)
  
  sex <- rep(NA_integer_, N)
  sex[ father ] <- 1
  sex[ mother ] <- 2

  L$father <- father
  L$mother <- mother
  L$sex <- sex
  L
}

