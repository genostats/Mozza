#' Kinship Matrix
#'
#' @param x a list of 'zygotes'
#'
#' @details If `x` is a list n zygotes, this function computes a kinship matrix 
#' based on the IBD sharing between haplotypes. The coefficients of the matrix
#' are always non-negative, and the diagonal coefficients are always greater than
#' or equal to one.

#' @return A \eqn{n \times n}{n x n} matrix 
#' @export
#'
#' @examples x1 <- zygote(100)
#' x2 <- zygote(100)
#' kinship.matrix( c(x1, x2) )
kinship.matrix <- function(x) kinship_matrix_(x)