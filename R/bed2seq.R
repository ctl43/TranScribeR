#' @export
constract_ideal_seq <- function(bed, ref_genome = Hsapiens, out){
  bed <- import(bed)
  bed <- .block2gr(bed)
  bed <- bed[seqnames(bed)%in%names(ref_genome)]
  p_bed <- bed[strand(bed)=="+"]
  p <- data.table(name=p_bed$name,seq=as.character(getSeq(ref_genome, p_bed)))
  p <- p[,list(paste0(seq, collapse = "")), by = name]
  m_bed <- bed[strand(bed)=="-"]
  m <- data.table(name=m_bed$name,seq=as.character(getSeq(ref_genome, m_bed)))
  m <- m[,list(paste0(rev(seq), collapse = "")), by = name]
  combined <- rbind(p, m)
  new_seq <- DNAStringSet(combined[[2]])
  names(new_seq) <- gsub(";.*","",combined[[1]])
  if(!is.null(out)){
    write_fasta(new_seq, out)
  }else{
    return(new_seq)
  }
}

.block2gr <- function(bed){
  new_start <- start(bed$blocks) + start(bed) - 1
  new_end <- end(bed$blocks) + start(bed) - 1
  new_names <- rep(bed$name, lengths(new_start))
  new_strand <- rep(strand(bed), lengths(new_start))
  new_seqnames <- rep(seqnames(bed), lengths(new_start))
  GRanges(new_seqnames, IRanges(unlist(new_start), unlist(new_end)), strand = new_strand, name=new_names)
}