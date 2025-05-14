#' Zygotes
#' 
#' @name zygote
#' @aliases zygote_prob
#' 
#' @usage zygote(ntiles, mean_length_tiles = 20)
#' @usage zygote_probs(probaTiles, mean_length_tiles = 20)
#' 
#' @param ntiles number of different tiles
#' @param probaTiles vector of probabilities for each tiles
#' @param mean_length_tiles mean length of the tiles in cM
#' 
#' @details \code{haplotype} generates a mosaic haplotype with tiles numbered from \code{0} to \code{(ntiles - 1)}. 
#' The haplotype number is drawn uniformly. \code{haplotype_probs} generates a mosaic haplotype
#' with tiles numbered from \code{0} to \code{length(probaTiles) - 1}, drawn according to the probabilities
#' in \code{probaTiles}. In both cases, the length of tiles is drawn in a exponential distribution with 
#' mean length \code{mean_length_tiles}.
#' 
#' @return an external pointer to an haplotype
#' @seealso \link{haplotype}, \link{zygote_peek}
#' @export zygote
#' @export zygote_probs
#' 
NULL

#' Print method for zygotes
#' @name print.zygote
#' @usage print(x, ...)
#' @exportS3Method print zygote
print.zygote <- function(x, ...) {
  if(is.list(x)) 
    cat("Vector of", length(x), "zygotes\n")
  else
    cat("Zygote")
}

#' Str method for zygotes
#' @name str.zygote
#' @usage str(object, ...)
#' @exportS3Method str zygote
str.zygote <- function(x, ...) {
  if(is.list(x)) 
    cat(" zyg [1:", length(x), "] zygotes\n", sep = "")
  else
    cat(" zygote")
}

