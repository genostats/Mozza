#' Modify a tile number
#'
#' @param x a single zygote 
#' @param tile.index.a a vector of index of tiles to modify on first haplotype
#' @param tile.a a vector of new tile numbers
#' @param tile.index.b a vector of index of tiles to modify on second haplotype
#' @param tile.b a vector of new tile numbers
#'
#' @details modify (in place) tiles in a zygote. Tiles are R-indexed.
#'
#' @return NULL
#' @export
zygote.poke <- function(x, tile.index.a, tile.a, tile.index.b, tile.b) {
  if(!isa(x, "zygote")) stop("Not an zygote")
  if(length(x) != 1) stop("Not a single zygote")
  zygote_poke(x[[1]], tile.index.a, tile.a, tile.index.b, tile.b)
}
