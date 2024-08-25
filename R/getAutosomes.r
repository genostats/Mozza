#' Size of autosomes
#'
#' @name getAutosomes
#' @aliases setAutosomes
#' 
#' @usage getAutosomes()
#' @usage setAutosomes(x)
#'
#' @param x a vector of lengthes (in centiMorgan)
#' 
#' @details getAutosomes() sends back the length of autosomes used by Mozza. At startup
#' it is a vector of length 22 corresponding to the human genome. setAutosomes() enables
#' to change this.
#' 
#' @return getAutosomes() returns a vector.
#' 
#' @examples # get current value
#' auto <- getAutosomes()
#' auto
#' # assuming now 2 chromosomes with length 100 and 10 cM
#' setAutosomes(c(100, 10))
#' # testing this on one individual
#' z <- zygote(50)
#' zygote_peek(z)
#' # setting back to the previous value
#' setAutosomes(auto)
#'
NULL
