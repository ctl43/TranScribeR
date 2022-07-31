#' @export
genome2tx <- function(target = NULL, identifer = "name", anno){
  names(anno) <- mcols(anno)[[identifer]]
  if(is.null(target)){
    site_1  <- anno$thickStart
    site_2  <- anno$thickEnd
  }else{
    anno <- anno[mcols(target)[identifer]]
    if(!all(strand(anno)==strand(target))){
      stop("The strand of target and annotation should be the same.")
    }
    site_1  <- start(target)
    site_2  <- end(target)
  }
  is_empty <- site_1==0|site_2==0
  # Gathering data
  strand <- as.character(strand(anno))
  is_m <- strand == "-"
  exon_len <- mcols(anno)$blockSizes
  exon_cs <- cumsum(exon_len)
  exon_start <- mcols(anno)$blockStarts + start(anno) - 1L
  exon_end <- exon_start + exon_len - 1L

  # Checking which exon the features falling at
  site_1_at <- unname(IRanges::which(site_1 >= exon_start & site_1 <= exon_end))
  site_2_at <- unname(IRanges::which(site_2 >= exon_start & site_2 <= exon_end))

  if(any(lengths(site_1_at)>1) | any(lengths(site_2_at)>1)){
    stop("More than one start codon/stop codon for each transcript")
  }

  # Computing the transcript coordination
  ee2site1_dist <- exon_end[site_1_at] - site_1
  es2site2_dist <- site_2 - exon_start[site_2_at]
  tx_site_1 <- exon_cs[site_1_at] - ee2site1_dist
  tx_site_2 <- exon_cs[site_2_at] - (exon_len[site_2_at] - es2site2_dist) + 1
  tx_start <- tx_site_1
  tx_start[is_m] <- sum(exon_len[is_m]) - tx_site_2[is_m]
  tx_stop <- tx_site_2
  tx_stop[is_m] <- sum(exon_len[is_m]) - tx_site_1[is_m] + 1
  tx_start[is_empty] <- 0
  tx_stop[is_empty] <- 0
  return(GRanges(seqnames = names(anno), IRanges(unlist(tx_start), unlist(tx_stop))))
}
