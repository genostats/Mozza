#' Mosaic Haplotypes
#' 
#' @name haplotype
#' @aliases haplotype.probs
#' 
#' @usage haplotype(ntiles, mean.length.tiles = 20)
#' @usage haplotype.probs(probaTiles, mean.length.tiles = 20)
#' 
#' @param ntiles number of different tiles
#' @param probaTiles vector of probabilities for each tiles
#' @param mean.length.tiles mean length of the tiles in cM
#' 
#' @details \code{haplotype} generates a mosaic haplotype with tiles numbered from \code{0} to \code{(ntiles - 1)}. 
#' The haplotype number is drawn uniformly. \code{haplotype.probs} generates a mosaic haplotype
#' with tiles numbered from \code{0} to \code{length(probaTiles) - 1}, drawn according to the probabilities
#' in \code{probaTiles}. In both cases, the length of tiles is drawn in a exponential distribution with 
#' mean length \code{mean.length.tiles}.
#' 
#' @seealso \link{zygote}, \link{haplotype.peek}
#' @return an external pointer to an haplotype
#' @export haplotype
#' @export haplotype.probs
#'

haplotype <- function(ntiles, mean.length.tiles = 20) {
  L <- list(haplotype_(ntiles, mean.length.tiles))
  class(L) <- "haplotype"
  L
}

haplotype.probs <- function(proba.tiles, mean.length.tiles = 20) {
  L <- list(haplotype_probs(proba.tiles, mean.length.tiles))
  class(L) <- "haplotype"
  L
}

#' Print method for haplotypes
#' @name print.haplotype
#' @usage print(x, ...)
#' @exportS3Method print haplotype
print.haplotype <- function(x, ...) {
  cat("Vector of", length(x), "haplotypes\n")
}

#' Str method for haplotypes
#' @name str.haplotype
#' @usage str(object, ...)
#' @exportS3Method str haplotype
str.haplotype <- function(x, ...) {
  cat(" hap [1:", length(x), "] haplotypes\n", sep = "")
}

# pour que le subsetting de vecteur d'haplotype reste un vecteur d'haplotypes
#' @exportS3Method `[` haplotype
`[.haplotype` <- function (x, ..., drop = FALSE) {
  y <- NextMethod("[")
  class(y) <- "haplotype"
  y
}

# Même problème avec la concaténation, qui ne propage pas la classe
#' @exportS3Method c haplotype
c.haplotype <- function (..., recursive = TRUE) {
  x <- list(...)
  if(any(sapply(x, \(z) !isa(z, "haplotype")))) 
    stop("Can't concatenate haplotypes with something else")
  y <- unlist(x, recursive = recursive)
  class(y) <- "haplotype"
  y
}


