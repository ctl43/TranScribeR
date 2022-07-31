#' @export
get_5putr <- function(cds, min_len = 50){
  fivep_utr_length <- start(cds) - 1
  pass <- fivep_utr_length >= min_len
  fivep_utr <- GRanges(seqnames(cds)[pass], IRanges(1, fivep_utr_length[pass]))
  fivep_utr$type <- "5prime_UTR"
  fivep_utr
}

#' @export
get_3putr <- function(cds, transcript_anno, min_len = 50, transcript_length = NULL){
  if(is.null(transcript_length)){
    len <- sum(transcript_anno[as.character(seqnames(cds))]$blockSizes)
  }else{
    len <- transcript_length
  }
  threep_utr_length <- len - end(cds)
  pass <- threep_utr_length >= min_len
  threep_utr <- GRanges(seqnames(cds)[pass], IRanges((len - threep_utr_length + 1)[pass], len[pass]))
  threep_utr$type <- "3prime_UTR"
  threep_utr
}
