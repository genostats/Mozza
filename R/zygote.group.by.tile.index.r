#' Group positions sharing the same tile index 
#'
#' @param x a single zygote 
#' @param chr vector of chromomose numbers
#' @param pos vector of position (in cM)
#'
#' @details sends two list giving the indexes of positions sharing the same tiles, for the two haplotypes
#'
#' @return a list
#' @export
zygote.group.by.tile.index <- function(x, chr, pos) {
  if(!isa(x, "zygote")) stop("Not an zygote")
  if(length(x) != 1) stop("Not a single zygote")
  L <- zygote_group_by_tile_index(x[[1]], chr, pos)
  names(L) <- c("a", "b")
  L
}
