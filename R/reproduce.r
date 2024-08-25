#' Mating
#'
#' @name reproduce
#'
#' @param z1,z2 zygotes
#'
#' @details If `z1` and `z2` are zyogtes, this function will proceed to
#' meiosis and send back an offspring.
#'
#' @return A zygote
#' @export
#'
#' @examples z1 <- zygote(100)
#' z2 <- zygote(100)
#' z3 <- reproduce(z1, z2)
#' kinship.matrix( c(z1, z2, z3) )
#'
NULL
