#' Generating unrelated individuals with tile haplotypes
#'
#' @param nb.inds number of (unrelated) individuals
#' @param haplos haplotype bed matrix
#' @param proba.haplos matrix of haplotypes probabilities
#' @param tile.length tile length in cM
#' @param kinship Logical. TRUE to get kinship matrix computed from IBD sharing.
#' @param fraternity Logical. TRUE to get fraternity matrix computed from IBD sharing.
#' 
#' @details If `proba.haplos` is missing, all haplotypes are used with the same probability.
#' `proba.haplos` should be a matrix with as many rows than there are haplotypes in `haplos`. 
#' If there are several columns in `proba.haplos`, 
#' `nb.inds` individuals will be generated for each column of probabilities. If `nb.inds` is 
#' a vector of length `ncol(proba.haplos)`, it specifies the number of individuals to be 
#' generated for each probability set.
#'
#' @return a list width components `bed`, `kinship` and `fraternity` (if applicable).
#' @export
#'
#' @examples #' # installs KGH if not already installed
#' if(!require("KGH")) install.packages("KGH", repos="https://genostats.github.io/R/")
#' # loads a bed matrix of 1006 european haplotypes
#' filepath <- system.file("extdata", "1KG_haplos.bed", package = "KGH")
#' H <- read.bed.matrix(filepath)
#' # Generates a bed matrix with 1000 individuals
#' x <- make.inds(1000, H)

make.inds <- function(nb.inds, haplos, proba.haplos, tile.length = 20, kinship = FALSE, fraternity = FALSE) {

  if(all(haplos@snps$dist == 0))
    stop("Set genetic distance between markers with set.dist !")

  if(missing(proba.haplos))
    L <- make_inds(nb.inds, tile.length, haplos@bed, haplos@snps$chr, haplos@snps$dist, kinship, fraternity)
  else {
    if(is.vector(proba.haplos))
      proba.haplos <- matrix(proba.haplos, ncol = 1)
    nb.inds <- rep_len(nb.inds, ncol(proba.haplos))
    L <- make_inds_probs(nb.inds, proba.haplos, tile.length, haplos@bed, haplos@snps$chr, haplos@snps$dist, kinship, fraternity)
  }

  ped <- data.frame(famid = 1:sum(nb.inds), id = 1:sum(nb.inds), father = NA, 
                    mother = NA, sex = NA, pheno = NA, stringsAsFactors = FALSE)
  
  x <- new("bed.matrix", bed = L$bed, snps = haplos@snps, ped = ped, p = NULL, mu = NULL,
           sigma = NULL, standardize_p = FALSE, standardize_mu_sigma = FALSE )
  if(getOption("gaston.auto.set.stats", TRUE))
    x <- set.stats(x)

  L$bed <- x
  L
}

