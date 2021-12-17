#' Generating unrelated individuals with tile haplotypes
#'
#' @param nb.inds number of (unrelated) individuals
#' @param haplos haplotype bed matrix
#' @param submap indices of SNPs with an effect
#' @param beta effect
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
#' set.seed(1)
#' I <- sample.int( ncol(H), 10000 )
#' beta <- rnorm(length(I))/sqrt(length(I))
#' x <- selected.inds(1000, H, I, beta)
#' X1 <- as.matrix( x$bed[,I] )
#' range( X1 %*% beta )

selected.inds <- function(nb.inds, haplos, submap, beta, proba.haplos, tile.length = 20, h2 = 0.5, kinship = FALSE, fraternity = FALSE) {

  if(all(haplos@snps$dist == 0))
    stop("Set genetic distance between markers with set.dist !")

  if(length(submap) != length(beta))
    stop("'submap' and 'beta' should have same length")

  if(any(submap < 1) | any(submap > ncol(haplos)))
    stop("bad submap")

  submap <- as.integer(submap) - 1L
  o <- order(submap)
  submap <- submap[o]
  beta   <- beta[o]

  if(missing(proba.haplos))
  L <- makeSelectedInds(nb.inds, tile.length, haplos@bed, haplos@snps$chr, haplos@snps$dist, submap, beta, h2, kinship, fraternity)
  else {
    stop("NOT YET IMPLEMENTED")
    if(is.vector(proba.haplos))
      proba.haplos <- matrix(proba.haplos, ncol = 1)
    nb.inds <- rep_len(nb.inds, ncol(proba.haplos))
    # L <- make_inds_probs(nb.inds, proba.haplos, tile.length, haplos@bed, haplos@snps$chr, haplos@snps$dist, kinship, fraternity)
  }

  # number of selected individuals
  nb.inds <- L$n
  L$n <- NULL
  ped <- data.frame(famid = 1:sum(nb.inds), id = 1:sum(nb.inds), father = NA, 
                    mother = NA, sex = NA, pheno = NA, stringsAsFactors = FALSE)
  
  x <- new("bed.matrix", bed = L$bed, snps = haplos@snps, ped = ped, p = NULL, mu = NULL,
           sigma = NULL, standardize_p = FALSE, standardize_mu_sigma = FALSE )
  if(getOption("gaston.auto.set.stats", TRUE))
    x <- set.stats(x)

  L$bed <- x
  L
}

