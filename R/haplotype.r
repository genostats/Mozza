#' Mosaic Haplotypes
#' 
#' @name haplotype
#' @aliases haplotype_probs
#' 
#' @usage haplotype(ntiles, mean_length_tiles = 20)
#' @usage haplotype_probs(probaTiles, mean_length_tiles = 20)
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
#' @seealso \link{zygote}, \link{haplotype_peek}
#' @return an external pointer to an haplotype
#' @export haplotype
#' @export haplotype_probs
#' 
NULL
