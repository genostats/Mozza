#' Kinship Matrix
#'
#' @param x a list of 'zygotes'
#' @param haplos a bed.matrix of haplotypes
#' @param filename file to be created
#'
#' @details If `x` is a list zygotes, creates a VCF file (with phased haplotypes)
#' 
#' @return NULL object
#' @export
#'
#' @examples #' # installs KGH if not already installed
#' if(!require("KGH")) install.packages("KGH", repos="https://genostats.github.io/R/")
#' # loads a bed matrix of 1006 european haplotypes
#' filepath <- system.file("extdata", "1KG_haplos.bed", package = "KGH")
#' H <- read.bed.matrix(filepath)
#' 
#' x <- make.inds(100, 1006)
#' filename <- tempfile(fileext = ".vcf")
#' drop.genotypes.vcf(x, H, filename)
drop.genotypes.vcf <- function(x, haplos, filename) {
  if(file.exists(filename)) stop("file already exists")
  # create VCF header
  z <- file(filename, "w")
  cat("##fileformat=VCFv4.1\n", file = z)
  cat('##FILTER=<ID=PASS,Description="All filters passed">\n', file = z)
  cat('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n', file = z)
  cat("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", file = z)
  n <- length(x)
  cat(paste0(sprintf("\tIND%0*d", floor(log10(n)) + 1, 1:n), collapse = ""), "\n", file = z)
  close(z)

  # append genotypes to VCF file
  drop_genotypes_to_vcf(x, filename, haplos@bed, haplos@snps$chr, haplos@snps$dist, 
                        haplos@snps$id, haplos@snps$pos, haplos@snps$A1, haplos@snps$A2)

  invisible(NULL)
}

