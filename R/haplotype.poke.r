#' Modify a tile number
#'
#' @param x a single haplotype 
#' @param tile.index a vector of indexes of tiles to modify
#' @param tile new tile numbers
#'
#' @details modify (in place) a tile in an haplotype. Tiles are R-indexed (started at 1).
#'
#' @return NULL
#' @export
haplotype.poke <- function(x, tile.index, tile) {
  if(!isa(x, "haplotype")) stop("Not an haplotype")
  if(length(x) != 1) stop("Not a single haplotype")
  haplotype_poke(x[[1]], tile.index, tile)
}
