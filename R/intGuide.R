#' @rdname intGuide
#' @title intGuide
#'
#' @author Cai Haodong
#'
#' @description Fetch on- and off-target from raw fastq files.
#'    This package integrated mergeReads, trimReads, alignBowtie2,
#'    sam2bam, bam2bed and guideBed, is used for a complete processing
#'    from raw fastq to sgRNA targets. Alternatively, you can also
#'    analysis your data step by step,in such mode you can tune
#'    lots of parameters to achieve better performance.
#'
#' @param input a folder containing raw fastq files. Not
#'     supported file path beginning with dot, as './'!!
#' @param ref genome reference for bowtie2
#' @param sgRNA sgRNA sequences with PAM
#' @param bowtie2 bowtie2 path
#' @param samtools samtools path
#' @param bedtools bedtools path
#'
#' @importFrom parallel detectCores
#'
#' @export
#'

intGuide <- function(input,
                     ref,
                     sgRNA,
                     bowtie2 = 'bowtie2',
                     samtools = 'samtools',
                     bedtools = 'bedtools'){

  # check sgRNA
  stopifnot('Provide sgRNA sequences with PAM!' = endsWith(sgRNA,'G')| endsWith(sgRNA,'g'))

  setwd(dirname(input))

  myfiles = dir(paste0('./',basename(input)) ,'q.gz$',full.names = TRUE)

  fileNameCore = unique(unlist(strsplit(basename(myfiles),'\\_R[1-2]'))[rep(c(TRUE,FALSE),length(myfiles)/2)])

  num_threads = parallel::detectCores()

  # raw data QC and merge reads
  if(!dir.exists('./1fastp')) dir.create('./1fastp')

  fq1_out = paste0('./1fastp/',fileNameCore,'_R1.fq')
  fq2_out = paste0('./1fastp/',fileNameCore,'_R2.fq')
  mg_out = paste0('./1fastp/',fileNameCore,'_merge.fq')

  for (i in 1:length(fileNameCore)) {
    mergeReads(fq1 = myfiles[2*i-1],
               fq2 = myfiles[2*i],
               fq1_out = fq1_out[i],
               fq2_out = fq2_out[i],
               mg_out = mg_out[i],
               threads = num_threads)
  }


  # trim reads
  if(!dir.exists('./2fasta')) dir.create('./2fasta')

  for (i in 1:length(fileNameCore)) {
    trimReads(fq1_out[i],fq2_out[i],mg_out[i],outdir = './2fasta')
  }


  # align
  if(!dir.exists('./3sam')) dir.create('./3sam')

  for (i in 1:length(fileNameCore)) {
    alignBowtie2(fa1 = paste0('./2fasta/',strsplit(basename(fq1_out[i]), "\\.")[[1]][1],'.fa'),
                 fa2 = paste0('./2fasta/',strsplit(basename(fq2_out[i]), "\\.")[[1]][1],'.fa'),
                 outdir = './3sam',
                 bowtie2 = 'bowtie2',
                 ref = ref,
                 threads = num_threads)
    alignBowtie2(fa1 = paste0('./2fasta/',strsplit(basename(mg_out[i]), "\\.")[[1]][1],'.fa'),
                 outdir = './3sam',
                 bowtie2 = 'bowtie2',
                 ref = ref,
                 threads = num_threads)
  }


  # sam to bam
  if(!dir.exists('./4bam')) dir.create('./4bam')

  for (i in 1:length(fileNameCore)) {
    sam2bam(sam = paste0('./3sam/',strsplit(basename(fq1_out[i]), "\\.")[[1]][1],'.sam'),
            outdir = './4bam',
            samtools = 'samtools')

    sam2bam(sam = paste0('./3sam/',strsplit(basename(mg_out[i]), "\\.")[[1]][1],'.sam'),
            outdir = './4bam',
            samtools = 'samtools')
  }

  # bam to bed

  if(!dir.exists('./5bed')) dir.create('./5bed')

  for (i in 1:length(fileNameCore)) {
    bam2bed(bam = paste0('./4bam/',strsplit(basename(fq1_out[i]), "\\.")[[1]][1],'.bam'),
            outdir = './5bed',
            bedtools = 'bedtools')

    bam2bed(bam = paste0('./4bam/',strsplit(basename(mg_out[i]), "\\.")[[1]][1],'.bam'),
            outdir = './5bed',
            bedtools = 'bedtools')
  }

  # bed to final res
  if(!dir.exists('./6guide')) dir.create('./6guide')

  for (i in 1:length(fileNameCore)) {
    guideBed(fqBed = paste0('./5bed/',strsplit(basename(fq1_out[i]), "\\.")[[1]][1],'.bed'),
             mgBed = paste0('./5bed/',strsplit(basename(mg_out[i]), "\\.")[[1]][1],'.bed'),
             outdir = './6guide',
             sgRNA = sgRNA)
  }

}
