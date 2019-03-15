
test_make_pairs <- function(pairs, len1, len2, SNPS = 1e4, haps = 100, tile.length = 20) {
  # il faut construire une bed matrix avec haps haplotypes différents
  Haplos <- as.bed.matrix( matrix( rbinom(haps*SNPS, 1, .5), nrow = haps) )
  chr <- sort( sample(1:22, SNPS, TRUE) )
  # positions au hasard
  dist <- as.vector(unlist(tapply(chr, chr, function(x) sort(runif(length(x), 0 , Human.autosomes.b37[x[1]])) ))) 
  # -- over --
  
  L <- make_pairs(pairs, len1, len2, tile.length, Haplos@bed, chr, dist)

  ids <- sprintf("m%0*d", log10(SNPS) + 1, 1:SNPS)
  snps <- data.frame(chr = chr, id = ids, dist = dist, pos = NA,
               A1 = NA, A2 = NA, stringsAsFactors = FALSE)

  ped <- data.frame(famid = rep(1:pairs, each=2), ids = c("A", "B"), father = NA, 
                    mother = NA, sex = NA, pheno = NA, stringsAsFactors = FALSE)

  x <- new("bed.matrix", bed = L$bed, snps = snps, ped = ped, p = NULL, mu = NULL,
      sigma = NULL, standardize_p = FALSE, standardize_mu_sigma = FALSE )

  list("bed.matrix" = x, "kin" = L$kin)
}

