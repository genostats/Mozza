#' Generating case control samples
#'
#' @param nb.inds number of (unrelated) individuals
#' @param haplos haplotype bed matrix
#' @param submap indices of SNPs with an effect
#' @param beta effect
#' @param proba.haplos matrix of haplotypes probabilities
#' @param proba.demes vector of demes probabilities / weights
#' @param tile.length tile length in cM
#' @param prevalence
#' @param h2
#' @param kinship Logical. TRUE to get kinship matrix computed from IBD sharing.
#' @param fraternity Logical. TRUE to get fraternity matrix computed from IBD sharing.
#' @param population vector giving the population of origin of each haplotype 
#' @param ... addditional parameters give the probability of drawing a tile in each population.
#' 
#' 
#' @details 
#' If provided, `proba.haplos` should be a matrix with as many rows than there are haplotypes in `haplos`. 
#' If there are several columns in `proba.haplos`, each column corresponds to a deme. The relative weights of the demes
#' in the total population is specified by `proba.demes`; if missing, all demes have same weight.
#' @details
#' If `proba.haplos` is missing, the `population` parameter and the additionnal 
#' parameters can be used to construct a probabiity matrix. Otherwise, all haplotypes are used with the same probability.
#' @details
#' The genotype effect `G.raw` is rescaled in `G.adj` to have variance `h2`. The environmental component `E` has variance `1 - h2`.
#' Both `G.adj` and `E` are centered unless `offset.demes` values are provided; in that case the values provided give the mean environmental
#' component in each deme, and `G.adj` is offsetted accordingly. The resulting liability `G.adj + E` is always centered with 
#' unit variance.
#'
#' @return a list width components `bed`, `kinship` and `fraternity` (if applicable).
#' @export
#'
#' @examples #' # installs KGH if not already installed
#' if(!require("KGH")) install.packages("KGH", repos="https://genostats.github.io/R/")
#' # loads a bed matrix of 1006 european haplotypes
#' filepath <- system.file("extdata", "1KG_haplos.bed", package = "KGH")
#' H <- read.bed.matrix(filepath)
#' 
#' # Set of causal SNPs / effects
#' set.seed(1)
#' I <- sample.int( ncol(H), 10000 )
#' beta <- rnorm(length(I))/sqrt(length(I))
#'
#' # first test. Prevalence 0, only controls : general population
#' # 1000 individuals in 11 demes with different proportions of TSI / IBS
#' p <- seq(0, 1, length = 11)
#' x.1 <- cases.controls(0, 1000, H, I, beta, prevalence = 0, TSI = p, IBS = 1 - p)
#' table(x.1$bed@ped$Deme)
#' z <- LD.thin(select.snps(x.1$bed, maf > 0.05), 0.1)
#' K <- GRM(z)
#' plot( eigen(K)$vectors, col = hsv(x.1$bed@ped$TSI) )
#'
#' # second test. Default prevalence and heritability, all haplotypes equal.
#' x.2 <- cases.controls(500, 500, H, I, beta)
#'
#' # third test
#' x.3 <- cases.controls(500, 500, H, I, beta, TSI = p, IBS = 1 - p)
#' 
#' # fourth test
#' x.4 <- cases.controls(500, 500, H, I, beta, offset.demes = seq(0,1,length=11), TSI = p, IBS = 1 - p)

cases.controls <- function(nb.cases, nb.controls, haplos, submap, beta, proba.haplos, proba.demes, offset.demes,
                                tile.length = 20, prevalence = 0.01, h2 = 0.5, kinship = FALSE, fraternity = FALSE,
                                population = haplos@ped$population, ...) {


  if(all(haplos@snps$dist == 0))
    stop("Set genetic distance between markers with set.dist !")

  if(length(submap) != length(beta))
    stop("'submap' and 'beta' should have same length")

  if(any(submap < 1) | any(submap > ncol(haplos)))
    stop("bad submap")

  probs <- list(...)

  if(missing(proba.haplos)) {
    proba.haplos <- NULL
  }
  if(length(probs) > 0 && !is.null(proba.haplos)) {
    warning("A matrix proba.haplos was given. Arguments ", names(probs), " are ignored.")
    probs <- list()
  }
  if(length(probs) > 0 && is.null(proba.haplos)) {
    proba.haplos <- make.proba.haplos(population, probs)
  }

  submap <- as.integer(submap) - 1L
  o <- order(submap)
  submap <- submap[o]
  beta   <- beta[o]

  threshold <- qnorm(1-prevalence)

  if(is.null(proba.haplos)) {
    L <- liabilitySelectedInds( c(nb.cases, nb.controls), c(threshold, -Inf), c(Inf, threshold), 
                           tile.length, haplos@bed, haplos@snps$chr, haplos@snps$dist, submap, beta, h2, kinship, fraternity)
  } else {
    if(is.vector(proba.haplos))
      proba.haplos <- matrix(proba.haplos, ncol = 1)
    if(missing(proba.demes)) 
      proba.demes <- rep(1/ncol(proba.haplos), ncol(proba.haplos))
    else
      proba.demes <- proba.demes / sum(proba.demes)
    if(missing(offset.demes))
      offset.demes <- rep(0, ncol(proba.haplos))

    L <- liabilitySelectedIndsProbs( c(nb.cases, nb.controls), proba.haplos, proba.demes, offset.demes, c(threshold, -Inf), c(Inf, threshold),
                                     tile.length, haplos@bed, haplos@snps$chr, haplos@snps$dist, submap, beta, h2, kinship, fraternity)
  }

  # number of selected individuals
  nb.inds <- L$n
  L$n <- NULL
  ped <- data.frame(famid = 1:sum(nb.inds), id = 1:sum(nb.inds), father = NA, 
                    mother = NA, sex = NA, pheno = c(rep(1,nb.cases), rep(0, nb.controls)), 
                    stringsAsFactors = FALSE)
  
  x <- new("bed.matrix", bed = L$bed, snps = haplos@snps, ped = ped, p = NULL, mu = NULL,
           sigma = NULL, standardize_p = FALSE, standardize_mu_sigma = FALSE )
  if(getOption("gaston.auto.set.stats", TRUE))
    x <- set.stats(x)

  x@ped$G.raw <- L$G.raw
  x@ped$G.adj <- L$G.adj
  x@ped$E <- L$E
  
  if(!is.null(proba.haplos)) {
    x@ped$Deme <- L$Deme + 1
  }
  if(length(probs) > 0) {
    for(a in names(probs)) {
      x@ped[[a]] <- probs[[a]][L$Deme + 1]
    }
  }

  L$bed <- x
  L
}

