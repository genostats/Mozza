#' Peek into haplotypes and zygotes
#' @name peek
#' @aliases haplotype_peek zygote_peek
#' 
#' @usage haplotype_peek(x)
#' @usage zygote_peek(x)
#' 
#' @param x a haplotype or a zygote
#' 
#' @details \code{haplotype_peek} sends back a data frame with columns 'chr', 'pos', 'tile'.
#' The position in 'pos' is the endpoint of the tile (in cM) whose number is in 'tile'. 
#' The last position is always the end of the chromosome.
#' \code{zygote_peek} sends back a list of two such data frames.
#' 
#' @examples 
#' # an example with long tiles
#' set.seed(1)
#' x <- haplotype(10, 100)
#' haplotype_peek(x)
NULL
