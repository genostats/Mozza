#' Mating
#'
#' @name reproduce
#'
#' @param z1,z2 zygote vectors of same length
#'
#' @details This function will create offsprings of individuals in z1 x z2
#'
#' @return A zygote vector
#' @export
#'
#' @examples
#' fa <- make.inds(2, 1000)
#' mo <- make.inds(2, 1000)
#' of <- reproduce(fa, mo)
#' K <- kinship.matrix( c(fa, mo, of) )
#' round(K, 1)
#' # You can also use `+` for a simple call. 
#' of2 <- fa + mo
#' L <- kinship.matrix( c(of, of2) )
#' round(L, 1)

reproduce <- function(z1, z2) {
  reproduce_vec(z1, z2)
}

#' @exportS3Method `+` zygote
`+.zygote` <- function(e1, e2) reproduce(e1, e2)

