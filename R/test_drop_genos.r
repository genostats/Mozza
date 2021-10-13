
test_drop_genos <- function(N = 10, chr.max = min(22,floor(N/3))) {
  Haplos <- as.bed.matrix( matrix( rbinom(4*N, 1, .5), nrow = 4) ) # une bed matrix avec haps haplotypes 
  chr <- sort( sample(1:chr.max, N, TRUE) )
  # positions au hasard
  dist <- as.vector(unlist(tapply(chr, chr, function(x) sort(runif(length(x), 0 , Human.autosomes.b37[x[1]])) ))) 

  bed <- test_xptr(Haplos@bed, chr, dist)

  id <- sprintf("m%0*d", log10(N) + 1, 1:N)
  snps <- data.frame(chr = chr, id = id, dist = dist, pos = NA,
               A1 = NA, A2 = NA, stringsAsFactors = FALSE)

  ped <- data.frame(famid = NA, id = c("F", "M", "O1", "O2"), father = c(NA,NA,"F","F"), 
                    mother = c(NA,NA,"M","M"), sex = NA, pheno = NA, stringsAsFactors = FALSE)

  new("bed.matrix", bed = bed, snps = snps, ped = ped, p = NULL, mu = NULL,
      sigma = NULL, standardize_p = FALSE, standardize_mu_sigma = FALSE )

}

