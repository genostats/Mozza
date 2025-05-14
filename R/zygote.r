#' Zygotes
#' 
#' @name zygote
#' @aliases zygote.prob
#' 
#' @usage zygote(ntiles, mean.length.tiles = 20)
#' @usage zygote.probs(probaTiles, mean.length.tiles = 20)
#' 
#' @param ntiles number of different tiles
#' @param proba.tiles vector of probabilities for each tiles
#' @param mean.length.tiles mean length of the tiles in cM
#' 
#' @details \code{zygote} generates a zygote of two mosaic haplotypes, with tiles numbered from \code{0} to \code{(ntiles - 1)}. 
#' The haplotype number is drawn uniformly. 
#' \code{zygote.probs} generates a zygote of two mosaic haplotypes
#' with tiles numbered from \code{0} to \code{length(probaTiles) - 1}, drawn according to the probabilities
#' in \code{probaTiles}. In both cases, the length of tiles is drawn in a exponential distribution with 
#' mean length \code{mean.length.tiles}.
#' 
#' @return an external pointer to an haplotype
#' @seealso \link{haplotype}, \link{zygote.peek}
#' @export zygote
#' @export zygote.probs
#' 
zygote <- function(ntiles, mean.length.tiles = 20) {
  L <- list(zygote_(ntiles, mean.length.tiles))
  class(L) <- "zygote"
  L
}

zygote.probs <- function(proba.tiles, mean.length.tiles = 20) {
  L <- list(zygote_probs(proba.tiles, mean.length.tiles))
  class(L) <- "zygote"
  L
}

#' Print method for zygotes
#' @name print.zygote
#' @usage print(x, ...)
#' @exportS3Method print zygote
print.zygote <- function(x, ...) {
  cat("Vector of", length(x), "zygotes\n")
}

#' Str method for zygotes
#' @name str.zygote
#' @usage str(object, ...)
#' @exportS3Method str zygote
str.zygote <- function(x, ...) {
  cat(" zyg [1:", length(x), "] zygotes\n", sep = "")
}

# pour que le subsetting de vecteur de zygote reste un vecteur de zygotes
#' @exportS3Method `[` zygote
`[.zygote` <- function (x, ..., drop = FALSE) {
  y <- NextMethod("[")
  class(y) <- "zygote"
  y
}

# Même problème avec la concaténation, qui ne propage pas la classe
#' @exportS3Method c zygote
c.zygote <- function (..., recursive = TRUE) {
  x <- list(...)
  if(any(sapply(x, \(z) !isa(z, "zygote")))) 
    stop("Can't concatenate zygotes with something else")
  y <- unlist(x, recursive = recursive)
  class(y) <- "zygote"
  y
}


