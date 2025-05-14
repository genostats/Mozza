#' Show haplotype composition
#'
#' @details sends a list of data frames, each with three column: `chr`, `pos`, and `tile`.
#' The position is in cM. Tile is the tile number (lowest possible value = 1).
#'
#' @return A list of data frames.
#' @export
haplotype.peek <- function(x) lapply(x, haplotype_peek)
