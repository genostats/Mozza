#' Group positions sharing the same tile index
#'
#' @param x a single haplotype 
#' @param chr vector of chromomose numbers
#' @param pos vector of position (in cM)
#'
#' @details sends a list giving the indexes of positions sharing the same tiles.
#'
#' @return a list
#' @export
haplotype.group.by.tile.index <- function(x, chr, pos) {
  if(!isa(x, "haplotype")) stop("Not an haplotype")
  if(length(x) != 1) stop("Not a single haplotype")
  haplo_group_by_tile_index(x[[1]], chr, pos)
}
