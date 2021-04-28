#' @title readBed
#'
#' @description Read bed file
#'
#' @param bed a .bed file
#'
#' @importFrom data.table fread
#'
#' @export

readBed <- function(bed){
  col.name = c('chr','start','end','name','quality','strand')
  col.class= c('character','numeric','numeric','character','numeric','character')
  fread(bed,col.names = col.name,colClasses = col.class,header = FALSE)
}


#' @title closeSeqIdx
#'
#' @description Find duplicates based on edit distance
#'
#' @param str_vec a character vector
#' @param len UMI length
#'
#' @export

closeSeqIdx <- function(str_vec,len){
  dup_idx = NULL
  for (i in 2:(len-1)) {
    for (j in (i+1):len) {
      dup_idx = c(dup_idx,which(duplicated(
        paste0(substr(str_vec,1,nchar(str_vec)-j),
               substr(str_vec,nchar(str_vec)-j+2,nchar(str_vec)-i),
               substr(str_vec,nchar(str_vec)-i+2,nchar(str_vec)))
      )))
    }
  }
  unique(dup_idx)
}


#' @title barcTrans
#'
#' @description Translate barcode sequences to number
#'
#' @param barc_vec a vector containing short sequences need to be parsed
#' @param barcode a vector containing the reference barcodes.
#'     Default is c('ATCACG','CGATGT','TTAGGC','TGACCA','ACAGTG',
#'     'GCCAAT','CAGATC','ACTTGA','GATCAG')
#'
#' @importFrom utils adist
#'
#' @export


barcTrans <- function(barc_vec,barcode = def_barcode){
  def_barcode = c('ATCACG','CGATGT','TTAGGC','TGACCA','ACAGTG',
                  'GCCAAT','CAGATC','ACTTGA','GATCAG')
  code = rep(0,length(barc_vec))
  for(i in 1:length(barcode)) code[which(adist(barc_vec,barcode[i])<2)] = i
  code
}


#' @title alignBowtie2
#'
#' @description Using bowtie2 to align reads
#'
#' @param fa1 a fasta file
#' @param fa2 a fasta file
#' @param outdir the folder of output
#' @param bowtie2 bowtie2 path
#' @param ref reference genome
#' @param threads number of threads used in alignment
#'
#' @export

alignBowtie2 <- function(fa1,
                         fa2 = NULL,
                         outdir,
                         bowtie2 = 'bowtie2',
                         ref,
                         threads = 8){
  # outdir check
  if(is.null(outdir)){
    outdir = dirname(fa1)
  }else stopifnot('Invalid outdir' = dir.exists(outdir))

  if(is.null(fa2)) {
    system(paste0(bowtie2,' -f --end-to-end -t --no-unal -p',threads,' -x ',ref,' -U ',fa1,' -S ',
                  paste0(outdir,'/',strsplit(basename(fa1),"\\.")[[1]][1],'.sam')))
  } else {
    system(paste0(bowtie2,' -f --end-to-end -t --no-unal --no-mixed --no-discordant --dovetail --no-contain --no-overlap  -p',
                  threads,' -x ',ref,' -1 ',fa1,' -2 ',fa2,' -S ',
                  paste0(outdir,'/',strsplit(basename(fa1),"\\.")[[1]][1],'.sam')))
  }
}


#' @title sam2bam
#'
#' @description Using samtools to transform sam to bam
#'
#' @param sam a .sam file
#' @param outdir the folder of output
#' @param samtools samtools path
#'
#' @export

sam2bam <- function(sam,
                    outdir,
                    samtools = 'samtools'){
  # outdir check
  if(is.null(outdir)){
    outdir = dirname(sam)
  }else stopifnot('Invalid outdir' = dir.exists(outdir))

  system(paste0(samtools,' view -Sb ',sam,'> ',
                paste0(outdir,'/',strsplit(basename(sam),"\\.")[[1]][1],'.bam')))
}


#' @title bam2bed
#'
#' @description Using bedtools to transform bam to bed
#'
#' @param bam a .bam file
#' @param outdir the folder of output
#' @param bedtools bedtools path
#'
#' @export

bam2bed <- function(bam,
                    outdir,
                    bedtools = 'bedtools'){
  # outdir check
  if(is.null(outdir)){
    outdir = dirname(bam)
  }else stopifnot('Invalid outdir' = dir.exists(outdir))

  system(paste0(bedtools,' bamtobed -i ',bam,'> ',
                paste0(outdir,'/',strsplit(basename(bam),"\\.")[[1]][1],'.bed')))
}



#' @title prepBed
#'
#' @description preprocess bed
#'
#' @param x data.frame from bed
#'
#' @importFrom dplyr mutate
#'
#' @export

prepBed <- function(x){
  x = filter(x,quality >= 20)

  if(all(endsWith(x$name,c('1','2')))){
    x = x[rep(c(TRUE,FALSE),nrow(x)/2)]
    x = mutate(x,name = substr(name,1,nchar(name)-2))
  }

  x = mutate(x,
             barcode0 = substr(name,nchar(name)-5,nchar(name)),
             barcode = barcTrans(barcode0),
             umi = substr(name,nchar(name)-17,nchar(name)-6),
             side = substr(name,nchar(name)-19,nchar(name)-19),
             insiteCode = paste0(chr,':',start,',',end,':',strand,':',barcode,':',umi))

  dupIdx = closeSeqIdx(x$insiteCode,12)

  x[-dupIdx]
}




