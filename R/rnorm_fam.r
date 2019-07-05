

rnorm.fam <- function(x, mean = 0, s2 = 1, sd = sqrt(s2)) {
  if(is(x, "bed.matrix"))
    x <- x@ped 
  if(!is(x, "data.frame"))
    stop("x should be a bed.matrix or a data.frame")
  id <- paste0(x$famid, ":", x$id)
  father <- paste0(x$famid, ":", x$father)
  mother <- paste0(x$famid, ":", x$mother)
  no.father <- !(father %in% id)
  no.mother <- !(mother %in% id)
  if(any(xor(no.father, no.mother)))
    stop("Familial structure issue")
  
  r <- rep(NA_real_, nrow(x))
  names(r) <- id
  r[ no.father ] <- rnorm( sum(no.father) )

  w <- is.na(r)
  while( any(w) ) {
    r[w] <- suppressWarnings( rnorm( sum(w), (r[father[w]] + r[mother[w]])/2, sqrt(.5)) )
    w <- is.na(r)
  }
  mean + sd*r 
}
