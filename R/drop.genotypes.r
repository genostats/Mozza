#' Kinship Matrix
#'
#' @param x a list of 'zygotes'
#' @param haplos a bed.matrix of haplotypes
#'
#' @details If `x` is a list zygotes, returns a bed.matrix 
#' constructed from the template haplotypes in `haplos`
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
#' Z <- vector("list", 100)
#' for(i in 1:100) Z[[i]] <- zygote(1006)
#' 
#' bm <- drop.genotypes(Z, H)
drop.genotypes <- function(x, haplos) {
  L <- drop_genotypes(x, haplos@bed, haplos@snps$chr, haplos@snps$dist)
  nb.inds <- length(x)
  ped <- data.frame(famid = 1:sum(nb.inds), id = 1:sum(nb.inds), father = NA,
                    mother = NA, sex = NA, pheno = NA, stringsAsFactors = FALSE)

  x <- new("bed.matrix", bed = L$bed, snps = haplos@snps, ped = ped, p = NULL, mu = NULL,
           sigma = NULL, standardize_p = FALSE, standardize_mu_sigma = FALSE )

  if(getOption("gaston.auto.set.stats", TRUE))
    x <- set.stats(x)

  x
}

