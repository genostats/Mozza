### pas sûr qu'il faille garder l'argument subset
vcf.gen.line <- function(genotypes, err.prob, dp1, dp2, subset) {

  if(!missing(subset)) genotypes <- genotypes[subset]

  L <- length(genotypes)

  if(dp1 != dp2) {
    mrd <- sample(dp1:dp2, L, replace = TRUE)
  } else {
    mrd <- rep(dp1, L)
  } 

  probs <- genotypes/2
  errs1 <- err.prob
  Cerrs1 <- (1 - errs1)
  
  ## Separate into two haplotypes without loss of generality, all hets have alternative allee on hap A1
  A1 <- probs; A1[A1 == 0.5] <- 1
  A2 <- probs; A2[A2 == 0.5] <- 0
  
  ## Probabilities of ref, alt err reads depending on the genotype.
  bA1p <- A1;   bA1p[A1 == 0] <- Cerrs1;       bA1p[A1 == 1] <- (1/3)*errs1
  bA2p <- A1;   bA2p[A1 == 0] <- (1/3)*errs1;  bA2p[A1 == 1] <- Cerrs1
  bA3p <- (2/3)*errs1
  
  BA1p <- A2;   BA1p[A2 == 0] <- Cerrs1;       BA1p[A2 == 1] <- (1/3)*errs1
  BA2p <- A2;   BA2p[A2 == 0] <- (1/3)*errs1;  BA2p[A2 == 1] <- Cerrs1
  BA3p <- (2/3)*errs1
  
  ## Number of ref, alt, and error reads
  rd  <- rpois(L, (0.5*mrd*bA1p)+(0.5*mrd*BA1p))
  rd2 <- rpois(L, (0.5*mrd*bA2p)+(0.5*mrd*BA2p))
  rd3 <- rpois(L, (0.5*mrd*bA3p)+(0.5*mrd*BA3p))
  NAs <- which(is.na(genotypes))
  rd[NAs] <- 0; rd2[NAs] <- 0; rd3[NAs] <- 0

  ## --- phred likelihoods
  ## First, "raw" likelihoods
  pg00 <- (Cerrs1)^rd * (errs1/3)^rd2 * (2*errs1/3)^rd3
  pg01 <- (0.5*(Cerrs1)+0.5*(errs1/3))^rd * (0.5*(Cerrs1)+0.5*(errs1/3))^rd2 * (2*errs1/3)^rd3
  pg11 <- (errs1/3)^rd * (Cerrs1)^rd2 * (2*errs1/3)^rd3
 
  ## ban zero values
  pg00 <- ifelse(pg00 == 0, .Machine$double.xmin, pg00)
  pg01 <- ifelse(pg01 == 0, .Machine$double.xmin, pg01)
  pg11 <- ifelse(pg11 == 0, .Machine$double.xmin, pg11)

  ## here are the phred scaled likelihood
  ph00 <- (-10)*log10(pg00)
  ph01 <- (-10)*log10(pg01)
  ph11 <- (-10)*log10(pg11)

  ## normalise a la GATK
  phtab <- cbind(ph00, ph01, ph11) - pmin(ph00, ph01, ph11)
  
  ## GQ is the diff between smallest and 2nd smallest, found with median after normalisation, bounded at 99 and always an integer
  GQ <- GQ(phtab)

  ## hard called genotypes from phtab [could be from ph00 pg01 pg1 as well...]
  genoN <- GT(phtab)
  genoN <- c("0/0", "0/1", "1/1")[genoN]
  
  ## Genotype Probabilities
  GP <- cbind(pg00, pg01, pg11)
  GP <- GP / rowSums(GP) # will be rounded to .001 later

  ## Dosages (rounded to 0.001)
  DS <- round( GP %*% 0:2, 3 )

  ## Now that dosages are computed, round GP to .001
  GP <- round(GP, 3)

  pasteVcfElts(genoN, rd, rd2, rd3, GQ, phtab, GP, DS)
}

