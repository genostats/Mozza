#' Generating unrelated individuals with tile haplotypes
#'
#' @param nb.inds number of (unrelated) individuals
#' @param ntiles number of distinct tiles
#' @param proba.haplos matrix of haplotypes probabilities
#' @param tile.length tile length in cM
#' 
#' @details If `proba.haplos` is missing, all haplotypes are used with the same probability.
#' `proba.haplos` should be a matrix with as many rows distinct tiles. 
#' If there are several columns in `proba.haplos`, 
#' `nb.inds` individuals will be generated for each column of probabilities. If `nb.inds` is 
#' a vector of length `ncol(proba.haplos)`, it specifies the number of individuals to be 
#' generated for each probability set.
#'
#' @return a zygote vector
#' @export
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
#' # to get a bed matrix, use drop.genotypes (here on SNPs 1 to 100)
#' xbed <- drop.genotypes(x, H[,1:100])

make.inds <- function(nb.inds, ntiles, proba.haplos, tile.length = 20) {

  if(missing(proba.haplos))
    x <- make_inds(nb.inds, ntiles, tile.length)
  else {
    if(is.vector(proba.haplos))
      proba.haplos <- matrix(proba.haplos, ncol = 1)
    nb.inds <- rep_len(nb.inds, ncol(proba.haplos))
    x <- make_inds_probs(nb.inds, proba.haplos, tile.length)
  }

  x
}

