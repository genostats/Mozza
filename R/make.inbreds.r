
make.inbreds <- function(nb.inbreds, haplos, f, a, hbd.length, non.hbd.length, tile.length = 20, segments = FALSE) {
  if(all(haplos@snps$dist == 0))
    stop("Set genetic distance between markers with set.dist !")
 
  if(missing(hbd.length) | missing(non.hbd.length)) {
    hbd.length <- 1/(a*(1-f))
    non.hbd.length <- 1/(a*f)
  } else {
    if(!missing(a) | !missing(f))
    stop("Specify either f and a, or hbd.length and non.hbd.length")
  }
   
  L <- make_inbreds(nb.inbreds, hbd.length, non.hbd.length, tile.length, 
                  haplos@bed, haplos@snps$chr, haplos@snps$dist, segments)

  ped <- data.frame(famid = rep(1:nb.inbreds, each=2), id = c("A", "B"), father = NA, 
                    mother = NA, sex = NA, pheno = NA, stringsAsFactors = FALSE)
  
  x <- new("bed.matrix", bed = L$bed, snps = haplos@snps, ped = ped, p = NULL, mu = NULL,
           sigma = NULL, standardize_p = FALSE, standardize_mu_sigma = FALSE )
  
  if(getOption("gaston.auto.set.stats", TRUE))
    x <- set.stats(x)
  
  R <- list("bed.matrix" = x, "inbred.coef" = L$inb)
  if(segments) R$segments <- L$segments
  R
}
