#' Generates mixed populuation
#'
#' @param nb.inds number of (unrelated) individuals per 'deme'
#' @param n.tiles number of distinct tiles
#' @param population vector of length n.tiles, giving the population of origin of each haplotype 
#' @param ... Additional parameters give the probability of drawing a tile in each population.
#' 
#' @details The number of 'demes' is given by the length of the vectors of probabilities.
#' There will be `nb.inds` individuals on each deme. `nb.inds` could be a vector giving 
#' different number of individuals for each deme.
#' 
#' @return a list width components `zygotes` and `probas` (giving for each zygote the 
#' probabilities that were used to draw the tiles)
#' @export
#'
#' @examples #' # installs KGH is not already installed
#' if(!require("KGH")) install.packages("KGH", repos="https://genostats.github.io/R/")
#' # loads a bed matrix of 1006 european haplotypes
#' filepath <- system.file("extdata", "1KG_haplos.bed", package = "KGH")
#' H <- read.bed.matrix(filepath)
#' 
#' # 100 individuals x 11 demes with different proportions of TSI / IBS
#' p <- seq(0, 1, length = 11)
#' x.1 <- mix.pop(100, n.tiles = nrow(H), population = H@ped$population, TSI = p, IBS = 1 - p)
#' # let's do a quick PCA
#' xbed.1 <- drop.genotypes(x.1$zygotes, H)
#' z.1 <- LD.thin(select.snps(xbed.1, maf > 0.05), 0.1)
#' K.1 <- GRM(z.1)
#' par(mfrow=c(1,2))
#' plot( eigen(K.1)$vectors, col = hsv(x.1$probas$TSI) )
#' 
#' # to generate a mixture of 4 populations on a 11 x 11 square
#' f <- function(x,y) c( (1-y)*c(x, 1-x), y*c(x, 1-x))
#' N <- 11
#' t <- rep(seq(0,1,length=N), N); u <- rep(seq(0,1,length=N), each = N)
#' pp <- t(mapply(f, t, u))
#' x.2 <- mix.pop(10, n.tiles = nrow(H), population = H@ped$population, TSI = pp[,1], 
#'                IBS = pp[,2], CEU = pp[,3], FIN = pp[,4])
#' # PCA
#' xbed.2 <- drop.genotypes(x.2$zygotes, H)
#' z.2 <- LD.thin(select.snps(xbed.2, maf > 0.05), 0.1)
#' K.2 <- GRM(z.2)
#' plot( eigen(K.2)$vectors, col = rgb(x.2$probas$TSI, x.2$probas$IBS, x.2$probas$CEU) )

mix.pop <- function(nb.inds, n.tiles, population, tile.length = 20, ...) {
  if(length(population) != n.tiles)
    stop("Dimensions mismatch")
  probs <- list(...)
  proba.haplos <- make.proba.haplos(population, probs)
  nb.inds <- rep_len(nb.inds, ncol(proba.haplos))
  x <- make.inds(nb.inds, n.tiles, proba.haplos, tile.length)
  PROBS <- list()
  for(a in names(probs)) {
    PROBS[[a]] <- as.vector(mapply(rep, probs[[a]], nb.inds))
  }
  list(zygotes = x, probas = PROBS)
}

