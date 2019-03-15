
test_families_of_4 <- function(fams = 10, SNPS = 100, haps = 20, chr.max = min(22,floor(SNPS/3))) {
  # il faut construire une bed matrix avec haps haplotypes différents
  Haplos <- as.bed.matrix( matrix( rbinom(haps*SNPS, 1, .5), nrow = haps) )
  chr <- sort( sample(1:chr.max, SNPS, TRUE) )
  # positions au hasard
  dist <- as.vector(unlist(tapply(chr, chr, function(x) sort(runif(length(x), 0 , Human.autosomes.b37[x[1]])) ))) 

  bed <- families_of_4(fams, Haplos@bed, chr, dist)

  ids <- sprintf("m%0*d", log10(SNPS) + 1, 1:SNPS)
  snps <- data.frame(chr = chr, id = ids, dist = dist, pos = NA,
               A1 = NA, A2 = NA, stringsAsFactors = FALSE)

  ped <- data.frame(famid = rep(1:fams, each=4), ids = c("F", "M", "O1", "O2"), father = c(NA,NA,"F","F"), 
                    mother = c(NA,NA,"M","M"), sex = NA, pheno = NA, stringsAsFactors = FALSE)

  new("bed.matrix", bed = bed, snps = snps, ped = ped, p = NULL, mu = NULL,
      sigma = NULL, standardize_p = FALSE, standardize_mu_sigma = FALSE )

}

