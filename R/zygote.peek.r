#' Show zygote composition
#'
#' @details sends a list of data frames, each with three column: `chr`, `pos`, and `tile`.
#' The position is in cM. Tile is the tile number (lowest possible value = 1).
#'
#' @return A list of data frames.
#' @export
zygote.peek <- function(x) {
  L <- lapply(x, zygote_peek)
  L <- unlist(L, recursive = FALSE)
  names(L) <- paste0( rep(1:length(x), each = 2), c("a", "b") )
  L
}
