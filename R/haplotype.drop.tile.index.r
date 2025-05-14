#' Get index of the tile at given positions
#'
#' @param x a single haplotype 
#' @param chr vector of chromomose numbers
#' @param pos vector of position (in cM)
#'
#' @details sends the tile indexes at the position given in chr, pos.
#'
#' @return an integer vector
#' @export
haplotype.drop.tile.index <- function(x, chr, pos) {
  if(!isa(x, "haplotype")) stop("Not an haplotype")
  if(length(x) != 1) stop("Not a single haplotype")
  haplo_drop_tile_index(x[[1]], chr, pos)
}
