#' Generating haplotypes
#'
#' @param nb.haps number of mosaic haplotypes to generate
#' @param haplos (template) haplotype bed matrix
#' @param proba.haplos matrix of haplotypes probabilities
#' @param tile.length tile length in cM
#' @param ibd Logical. TRUE to get a matrix giving the proportion genome shared IBD by two haplotypes.
#' 
#' @details If `proba.haplos` is missing, all haplotypes from `haplos` are used with the same probability.
#' `proba.haplos` should be a matrix with as many rows than there are haplotypes. If there are several columns, 
#' in `proba.haplos`, 
#' `nb.haps` new haplotypes will be generated for each column of probabilities. If `nb.haps` is 
#' a vector of length `ncol(proba.haplos)`, it specifies the number of haplotypes to be 
#' generated for each probability set.
#'
#' @return a list width components `bed`, and `ibd` (if applicable).
#' @export
#'
#' @examples #' # installs KGH if not already installed
#' if(!require("KGH")) install.packages("KGH", repos="https://genostats.github.io/R/")
#' # loads a bed matrix of 1006 european haplotypes
#' filepath <- system.file("extdata", "1KG_haplos.bed", package = "KGH")
#' H <- read.bed.matrix(filepath)
#' # Generates a bed matrix with 1000 haplotypes
#' x <- make.haps(1000, H)

make.haps <- function(nb.haps, haplos, proba.haplos, tile.length = 20, ibd = FALSE) {

  if(all(haplos@snps$dist == 0))
    stop("Set genetic distance between markers with set.dist !")

  if(missing(proba.haplos))
    L <- make_haps(nb.haps, tile.length, haplos@bed, haplos@snps$chr, haplos@snps$dist, ibd)
  else {
    if(is.vector(proba.haplos))
      proba.haplos <- matrix(proba.haplos, ncol = 1)
    nb.haps <- rep_len(nb.haps, ncol(proba.haplos))
    L <- make_haps_probs(nb.haps, proba.haplos, tile.length, haplos@bed, haplos@snps$chr, haplos@snps$dist, ibd)
  }

  ped <- data.frame(famid = 1:sum(nb.haps), id = 1:sum(nb.haps), father = NA, 
                    mother = NA, sex = NA, pheno = NA, stringsAsFactors = FALSE)
  
  x <- new("bed.matrix", bed = L$bed, snps = haplos@snps, ped = ped, p = NULL, mu = NULL,
           sigma = NULL, standardize_p = FALSE, standardize_mu_sigma = FALSE )
  if(getOption("gaston.auto.set.stats", TRUE))
    x <- set.stats(x)

  L$bed <- x
  L
}

