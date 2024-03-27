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

  zz <- file(filename, "w")
  cat("##fileformat=VCFv4.1\n##fileDate=2021-10-12\n##source=GastonVCFsimulation\n", file=zz)
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
    cat("\t999\tPASS\t.\tGT:AD:DP:GQ:PL\t", file = zz)
    cat(vcf.gen.line(geno, eps[i], depth1, depth2, subset), file = zz)
    cat("\n", file = zz)
    if( i %% 100 == 0 ) cat("Writing SNP", i, "\r")
  }
  cat("\n")
  close(zz)
}

