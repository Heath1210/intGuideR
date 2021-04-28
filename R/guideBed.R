#' @rdname guideBed
#' @title guideBed
#'
#' @description Fetch on- and off-targets of sgRNA from bed
#'
#' @param fqBed the file path to fq bed file.
#' @param mgBed the file path to mg bed file.
#' @param outdir the folder of output
#' @param sgRNA sgRNA sequences with PAM
#'
#' @importFrom dplyr `%>%`
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom dplyr arrange
#' @importFrom dplyr summarize
#' @importFrom dplyr desc
#' @importFrom dplyr group_by
#' @importFrom data.table fwrite
#' @importFrom data.table fread
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg38
#' @importFrom motifRG getSequence
#' @importFrom GenomicRanges findOverlaps
#' @importFrom Biostrings pairwiseAlignment
#' @importFrom dplyr n
#'
#' @export
#'


guideBed <- function(fqBed,
                     mgBed,
                     outdir = NULL,
                     sgRNA){
  # outdir check
  if(is.null(outdir)){
    outdir = paste0(dirname(mgBed),'/result')
    if(!dir.exists(outdir)) dir.create(outdir)
  }else stopifnot('Invalid outdir' = dir.exists(outdir))

  fqbed = prepBed(readBed(fqBed))
  mgbed = prepBed(readBed(mgBed))

  bed = rbind(fqbed,mgbed)

  bed = filter(bed,chr %in% c(1:22,'MT','X','Y'))

  bed_grange =  GRanges(
    seqnames = bed$chr,
    ranges = IRanges(start = ifelse(bed$strand=='+',bed$start-10,bed$end-10),
                     end = ifelse(bed$strand=='+',bed$start+10,bed$end+10)),
    strand = bed$strand)

  bed_grange2 = GRanges(
    seqnames = bed$chr,
    ranges = IRanges(start = ifelse(bed$strand=='+',bed$start-10,bed$end-10),
                     end = ifelse(bed$strand=='+',bed$start+10,bed$end+10)),
    strand = ifelse(bed$strand=='+','-','+'))

  bed_overlap = findOverlaps(bed_grange,bed_grange2)

  bed = bed[unique(bed_overlap@from)]

  cut_grange = getGRange(bed)

  sgRNA = DNAString(sgRNA)

  sgRNA_align_f = pairwiseAlignment(getSequence(cut_grange, BSgenome.Hsapiens.UCSC.hg38),
                                    sgRNA,type = 'local-global')

  sgRNA_align_r = pairwiseAlignment(getSequence(cut_grange, BSgenome.Hsapiens.UCSC.hg38),
                                    reverseComplement(sgRNA),type = 'local-global')

  sgRNA_match_idx = c(which(sgRNA_align_f@score > 0),which(sgRNA_align_r@score > 0))

  deviding_idx = length(which(sgRNA_align_f@score>0))

  bed = bed[sgRNA_match_idx]

  cut_grange = getGRange(bed)

  sgRNA_align_f_only = pairwiseAlignment(getSequence(cut_grange[1:deviding_idx], BSgenome.Hsapiens.UCSC.hg38),
                                          sgRNA,type = 'local-global')

  sgRNA_align_r_only = pairwiseAlignment(getSequence(cut_grange[(deviding_idx+1):length(cut_grange)], BSgenome.Hsapiens.UCSC.hg38),
                                          reverseComplement(sgRNA),
                                          type = 'local-global')

  preciseGR_f = preciseGRange(bed[1:deviding_idx],sgRNA_align_f_only@pattern@range)

  preciseGR_r = preciseGRange(bed[(deviding_idx+1):nrow(bed)],sgRNA_align_r_only@pattern@range)

  bed = mutate(bed,
    pattern = c(as.character(sgRNA_align_f_only@pattern),
                as.character(reverseComplement(DNAStringSet(sgRNA_align_r_only@pattern)))),
    coordinate = paste0(c(preciseGR_f@seqnames,preciseGR_r@seqnames),':',
                        c(preciseGR_f@ranges@start,preciseGR_r@ranges@start))
  )

  align_stat = summarize(group_by(bed,pattern,coordinate),n())
  colnames(align_stat) <- c('pattern','coordinate','number')

  fwrite(align_stat,paste0(outdir,'/',strsplit(basename(fqBed), "\\.")[[1]][1],'_stat.csv'))
  fwrite(bed,paste0(outdir,'/',strsplit(basename(fqBed), "\\.")[[1]][1],'_full.csv'))
}



#' @title getGRange
#'
#' @description create GRange object from bed
#'
#' @param x a data.frame
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#'

getGRange <- function(x){
  GRanges(seqnames =ifelse(x$chr=='MT','chrM',paste0('chr',x$chr)),
          ranges = IRanges(start = ifelse(x$strand=='+',x$start-50,x$end-50),
                           end = ifelse(x$strand=='+',x$start+50,x$end+50)),
          strand = x$strand)
  }


#' @title preciseGRange
#'
#' @description get precise ranges from relative ranges
#'
#' @param x a data.frame
#' @param alignRange  IRange object
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom BiocGenerics end
#'

preciseGRange = function(x,alignRange){
  GRanges(seqnames =ifelse(x$chr=='MT','chrM',paste0('chr',x$chr)),
          ranges = IRanges(start = ifelse(x$strand=='+',x$start+alignRange@start-51,x$end-end(alignRange)+51),
                           end = ifelse(x$strand=='+',x$start+end(alignRange)-51,x$end-alignRange@start+51)),
          strand = x$strand)
}
