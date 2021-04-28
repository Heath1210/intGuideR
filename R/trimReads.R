#' @rdname trimReads
#' @title Reads filtration and trimming
#'
#' @description Filter reads matching primer and linker, dump others.
#'    Adjust reads to same direction that linker is on the left side.
#'    Trim off linker and primer sequence followed by adding UMI, barcode
#'    to reads ID which has been simplied but informative enough to
#'    distinguish each reads. Save fasta files for further analysis.
#'
#' @param path2fq1 the file path to fq1 file.
#' @param path2fq2 the file path to fq2 file.
#' @param path2mg the file path to merged fastq file.
#' @param outdir the file fold where output files put in
#' @param dsODN1 the sequence of dsODN closest to integrated host genome with
#'     18bp length recommended (forward)
#' @param dsODN2 the sequence of dsODN closest to integrated host genome with
#'     18bp length recommended (reverse)
#' @param linker the sequence of adaptor linker at the end, near the genome part.
#'     Default linker is modified from INSPIIRED pipeline.
#'     If you use other linker, change to corresponding sequence.
#'
#' @importFrom Biostrings vcountPattern
#' @importFrom Biostrings vmatchPattern
#' @importFrom Biostrings reverseComplement
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings DNAString
#' @importFrom ShortRead readFastq
#' @importFrom ShortRead writeFasta
#' @importFrom Biostrings width
#' @importFrom dplyr `%>%`
#'
#' @export
#'

trimReads <- function(path2fq1,
                      path2fq2,
                      path2mg,
                      outdir = NULL,
                      dsODN1 = 'ATATGTTAATAACGGTAT',
                      dsODN2 = 'ATGACAACTCAATTAAAC',
                      linker = 'CTCCGCTTAAGGGACT'
                      ){
  # check validity of input arguments
  stopifnot('Error: input file not exists, please check your file path or file name!' =
              file.exists(path2fq1) & file.exists(path2fq2) & file.exists(path2mg))

  # output
  if(is.null(outdir)){
    outdir = paste0(dirname(path2fq1),'/fasta')
    if(!dir.exists(outdir)) dir.create(outdir)
  }

  dsODN1 = DNAString(dsODN1)
  dsODN2 = DNAString(dsODN2)
  linker = DNAString(linker)

  # loading fastq
  mg = readFastq(path2mg)
  fq1 = readFastq(path2fq1)
  fq2 = readFastq(path2fq2)

  # QC total reads
  total_reads = length(fq1)+length(mg)
  merge_percentage = length(mg)/total_reads

  # find out the merged reads that can be matched with reverse and complement dsODN sequence
  mg_rc_idx1 = which(vcountPattern(reverseComplement(dsODN1),substr(mg@sread,width(mg)-40,width(mg)))==1)
  mg_rc_idx2 = which(vcountPattern(reverseComplement(dsODN2),substr(mg@sread,width(mg)-40,width(mg)))==1)

  stopifnot('Please adjust your strictness of matching!' =
              length(intersect(mg_rc_idx1,mg_rc_idx2))==0)

  # replace reverse and complement reads so that all reads have same orientation
  mg@sread[c(mg_rc_idx1,mg_rc_idx2)] = reverseComplement(mg@sread[c(mg_rc_idx1,mg_rc_idx2)])

  # positively filter reads that can be matched with dsODN and linker
  mg_match_dsODN1 = which(vcountPattern(dsODN1,substr(mg@sread,1,40))==1)
  mg_match_dsODN2 = which(vcountPattern(dsODN2,substr(mg@sread,1,40))==1)

  mg_match_linker = which(vcountPattern(reverseComplement(linker),
                                        substr(mg@sread,width(mg)-60,width(mg)-20))==1)

  mg_match_idx1 = intersect(mg_match_dsODN1,mg_match_linker)
  mg_match_idx2 = intersect(mg_match_dsODN2,mg_match_linker)

  mg = mg[c(mg_match_idx1,mg_match_idx2)]

  # get the human genome sequence start location in each reads
  mg_start_dsODN1 = unlist(vmatchPattern(dsODN1,substr(mg@sread,1,40))@ends)+1

  mg_start_dsODN2 = unlist(vmatchPattern(dsODN2,substr(mg@sread,1,40))@ends)+1

  mg_start = c(mg_start_dsODN1,mg_start_dsODN2)

  # get the human genome sequence end location in each reads
  mg_end = unlist(vmatchPattern(reverseComplement(linker),
                                    substr(mg@sread,width(mg)-60,
                                           width(mg)-20))@ends)+width(mg)-length(linker)-61

  # the start and end locations are calculated in a vectorization way, normally have same length
  # if not, maybe dsODN or linker matched twice in one reads, so the match boundaries should be narrowed
  stopifnot('Please adjust your strictness of matching!'= length(mg_start)==length(mg_end))

  # trim off dsOND and linker sequences in reads, keep only human genome sequences
  mg_trim = substr(mg@sread,mg_start,mg_end)

  # name trimed reads the shortest unique strings of reads IDs appending with UMI sequencs
  # normally reads IDs have two patterns: with and without barcode
  # the length of reads ID with barcode is larger than 60
  mg_id = data.frame(matrix(unlist(strsplit(as.character(mg@id),' ')),ncol = 3,byrow = T),stringsAsFactors = FALSE)[,1]

  names(mg_trim) = paste0(substr(mg_id,24,width(mg_id)),':',
                           c(rep('A',length(mg_match_idx1)),rep('B',length(mg_match_idx2))),':',
                           substr(mg@sread,mg_end+length(linker)+1,mg_end+length(linker)+18) %>%
                             DNAStringSet() %>% reverseComplement())

  mg_trim = mg_trim[which(nchar(mg_trim)>30)]

  # find out the fq1 reads that can be matched with reverse and complement dsODN sequence
  fq_rc_idx1 = which(vcountPattern(dsODN1,substr(fq2@sread,1,40))==1)
  fq_rc_idx2 = which(vcountPattern(dsODN2,substr(fq2@sread,1,40))==1)

  stopifnot('Please adjust your strictness of matching!' =
              length(intersect(fq_rc_idx1,fq_rc_idx2))==0)


  # exchange reads so that all reads are in same orientation
  fq1_temp = fq1
  fq1@sread[c(fq_rc_idx1,fq_rc_idx2)] = fq2@sread[c(fq_rc_idx1,fq_rc_idx2)]
  fq1@quality@quality[c(fq_rc_idx1,fq_rc_idx2)] = fq2@quality@quality[c(fq_rc_idx1,fq_rc_idx2)]

  fq2@sread[c(fq_rc_idx1,fq_rc_idx2)] = fq1_temp@sread[c(fq_rc_idx1,fq_rc_idx2)]
  fq2@quality@quality[c(fq_rc_idx1,fq_rc_idx2)] = fq1_temp@quality@quality[c(fq_rc_idx1,fq_rc_idx2)]


  # the index of PE reads that can be matched with dsODN and linker
  fq_match_dsODN1 = which(vcountPattern(dsODN1,substr(fq1@sread,1,40))==1)
  fq_match_dsODN2 = which(vcountPattern(dsODN2,substr(fq1@sread,1,40))==1)
  fq_match_linker = which(vcountPattern(linker,substr(fq2@sread,20,60))==1)

  fq_match_idx1 = intersect(fq_match_dsODN1,fq_match_linker)
  fq_match_idx2 = intersect(fq_match_dsODN2,fq_match_linker)

  # filter
  fq1 = fq1[c(fq_match_idx1,fq_match_idx2)]
  fq2 = fq2[c(fq_match_idx1,fq_match_idx2)]

  reads_matching_primers = length(mg)+length(fq1)

  # get genome sequence start location in fq1 reads
  fq_start_dsODN1 = unlist(vmatchPattern(dsODN1,substr(fq1@sread,1,40))@ends)+1
  fq_start_dsODN2 = unlist(vmatchPattern(dsODN2,substr(fq1@sread,1,40))@ends)+1

  fq_start = c(fq_start_dsODN1,fq_start_dsODN2)

  # get the human genome sequence end location in fq2 reads
  fq_end = unlist(vmatchPattern(linker,substr(fq2@sread,20,60))@ends)+20

  # matching results control
  stopifnot('Please adjust your strictness of matching!' =
              length(fq_start)==length(fq_end))

  # trim off dsODN and linker sequences, get genome sequences and UMI
  fq1_trim = substr(fq1@sread,fq_start,width(fq1))
  fq2_trim = substr(fq2@sread,fq_end,width(fq2))

  # name trimmed reads the shortest unique strings of reads IDs appending with UMI sequences
  names(fq1_trim) = names(fq2_trim) =
    paste0(substr(fq1@id,24,width(fq1@id)-ifelse(width(fq1@id)>60,24,15)),':',
           c(rep('A',length(fq_match_idx1)),rep('B',length(fq_match_idx2))),':',
           substr(fq2@sread,fq_end-length(linker)-18,fq_end-length(linker)-1))

  short_idx = which(nchar(fq1_trim)<20 | nchar(fq2_trim)<20)

  if(length(short_idx)>0){
    fq1_trim = fq1_trim[-short_idx]
    fq2_trim = fq2_trim[-short_idx]
  }

  reads_usable_for_sites_detection = length(fq1_trim)+length(mg_trim)

  writeFasta(DNAStringSet(mg_trim),
             paste0(outdir,'/',strsplit(basename(path2mg), "\\.")[[1]][1],'.fa'))

  writeFasta(DNAStringSet(fq1_trim),
             paste0(outdir,'/',strsplit(basename(path2fq1), "\\.")[[1]][1],'.fa'))

  writeFasta(DNAStringSet(fq2_trim),
             paste0(outdir,'/',strsplit(basename(path2fq2), "\\.")[[1]][1],'.fa'))

  log_info = list(
    total_reads = total_reads,
    merge_percentage = merge_percentage,
    reads_matching_primers = reads_matching_primers,
    reads_usable_for_sites_detection = reads_usable_for_sites_detection)

  fwrite(log_info,
         paste0(outdir,'/',strsplit(basename(path2fq1), "\\.")[[1]][1],'.log'),
         sep = '\t' )

 }
