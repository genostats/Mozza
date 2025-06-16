#' Generates VCF file
#'
#' @param x 
#' @param filename 
#' @param depth1 
#' @param depth2 
#' @param eps 
#' @export
#'
#' @examples
#' x <- as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)
#' vcf.gen(x, "/tmp/LCT.vcf", 5, 10)

vcf.gen <- function(x, filename, depth1, depth2, eps = 10^((-0.1)*runif(ncol(x), 20, 30)), subset) {
  if(!missing(subset)) 
    if(any(subset < 0) | any(subset > nrow(x))) 
      stop("Bad subset") 

  if(missing(depth1) | missing(depth2)) stop("Specify depth1 and depth2")

  zz <- file(filename, "w")
  cat("##fileformat=VCFv4.1\n##fileDate=", format(Sys.time(), "%Y%M%d"), "\n##source=Mozza\n", sep = "", file = zz)
  cat('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n', file = zz)
  cat('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">\n', file = zz)
  cat('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n', file = zz)
  cat('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality. Phred-scaled confidence that the genotype assignment is correct">\n', file = zz)
  cat('##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes">\n', file = zz)
  cat('##FORMAT=<ID=GP,Number=G,Type=Float,Description="Genotype posterior probabilities">\n', file = zz)
  cat('##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage of the alternate allele">\n', file = zz)
  cat("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", file = zz)
  
  if(anyDuplicated(x@ped$id)) {
    id <- paste(x@ped$famid, x@ped$id, sep=":")
  } else {
    id <- x@ped$id
  }
  if(missing(subset)) 
    cat(paste(id, collapse = "\t"), "\n", sep = "", file = zz)
  else
    cat(paste(id[subset], collapse = "\t"), "\n", sep = "", file = zz)
  
  num.chr <- as.integer(x@snps$chr)
  CHR <- ifelse(is.na(num.chr), as.character(x@snps$chr), as.character(num.chr))
  eps <- rep_len(eps, ncol(x))

  for(i in 1:ncol(x)) {
    geno <- get.geno.vector(x, i);
    cat(CHR[i], x@snps$pos[i], x@snps$id[i], x@snps$A1[i], x@snps$A2[i], file = zz, sep = "\t")
    cat("\t999\tPASS\t.\tGT:AD:DP:GQ:PL:GP:DS\t", file = zz)
    cat(vcf.gen.line(geno, eps[i], depth1, depth2, subset), file = zz)
    cat("\n", file = zz)
    if( i %% 100 == 0 ) cat("Writing SNP", i, "\r")
  }
  cat("\n")
  close(zz)
}

