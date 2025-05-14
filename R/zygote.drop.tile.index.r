#' Get index of the tile at given positions
#'
#' @param x a single zygote 
#' @param chr vector of chromomose numbers
#' @param pos vector of position (in cM)
#'
#' @details sends the tile indexes at the position given in chr, pos.
#'
#' @return a list of two integer vectors, for the two haplotypes
#' @export
zygote.drop.tile.index <- function(x, chr, pos) {
  if(!isa(x, "zygote")) stop("Not an zygote")
  if(length(x) != 1) stop("Not a single zygote")
  L <- zygote_drop_tile_index(x[[1]], chr, pos)
  names(L) <- c("a", "b")
  L
}
