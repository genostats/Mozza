#' Kinship Matrix
#'
#' @param x a list of 'zygotes'
#' @param haplos a bed.matrix of haplotypes
#' @param phased a boolean 
#'
#' @details If `x` is a list zygotes, returns a bed.matrix 
#' constructed from the template haplotypes in `haplos`. If phased is TRUE,
#' instead of genotypes, two haplotypes are built for each individual. 
#' 
#' @return A bed.matrix
#' @export
#'
#' @examples #' # installs KGH if not already installed
#' if(!require("KGH")) install.packages("KGH", repos="https://genostats.github.io/R/")
#' # loads a bed matrix of 1006 european haplotypes
#' filepath <- system.file("extdata", "1KG_haplos.bed", package = "KGH")
#' H <- read.bed.matrix(filepath)
#' 
#' x <- make.inds(100, 1006)
#' 
#' bm <- drop.genotypes(x, H)
drop.genotypes <- function(x, haplos, phased = FALSE) {
  L <- drop_genotypes(x, haplos@bed, haplos@snps$chr, haplos@snps$dist, phased)
  nb.inds <- length(x)

  if(phased) {
    ped <- data.frame(famid = rep(1:nb.inds, each = 2), id = paste0(rep(1:nb.inds, each = 2), c("_1", "_2")),
                      father = NA, mother = NA, sex = NA, pheno = NA, stringAsFactors = FALSE)
  } 
  else {
    ped <- data.frame(famid = 1:nb.inds, id = 1:nb.inds, father = NA,
                      mother = NA, sex = NA, pheno = NA, stringsAsFactors = FALSE)
  }

  x <- new("bed.matrix", bed = L$bed, snps = haplos@snps, ped = ped, p = NULL, mu = NULL,
           sigma = NULL, standardize_p = FALSE, standardize_mu_sigma = FALSE )

  if(getOption("gaston.auto.set.stats", TRUE))
    x <- set.stats(x)

  x
}

