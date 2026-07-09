#' HBD segments
#'
#' @param x a list of 'zygotes'
#'
#' @details If `x` is a list n zygotes, this function compute for each individual
#' the list of HBD segments along its genome.

#' @return A list of data frames.
#' @export
#'
#' @examples set.seed(1)
#' x <- zygote(5, 100) # only few founder haplotypes, with long tiles
#' HBD.segments(x)
HBD.segments <- function(x) HBD_segments_(x)
