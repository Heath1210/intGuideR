#' @title mergeReads
#'
#' @description Using fastp to merge reads
#'
#' @param fq1 a fastq file
#' @param fq2 a fastq file
#' @param fastp fastp path
#' @param fq1_out filtered fq1 output file
#' @param fq2_out filtered fq2 output file
#' @param mg_out merged reads output file
#' @param threads number of threads used in alignment
#'
#' @export

mergeReads <- function(fq1,
                       fq2,
                       fastp = 'fastp',
                       fq1_out = NULL,
                       fq2_out = NULL,
                       mg_out = NULL,
                       threads = 8){
  if(is.null(fq1_out)){
    fq1_out = paste0(dirname(fq1),'/',strsplit(basename(fq1),"\\.")[[1]][1],'.fq')
  }
  if(is.null(fq2_out)){
    fq2_out = paste0(dirname(fq2),'/',strsplit(basename(fq2),"\\.")[[1]][1],'.fq')
  }
  if(is.null(mg_out)){
    mg_out = paste0(dirname(fq1),'/',strsplit(basename(fq1),"\\.")[[1]][1],'_merge.fq')
  }

  system(paste0(fastp," -g -w ",threads," -i ",fq1," -o ",fq1_out,
                " -I ",fq2," -O ",fq2_out," -m --merged_out ",mg_out))

}

