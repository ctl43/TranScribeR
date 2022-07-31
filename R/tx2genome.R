#' @export
tx2genome <- function(tx_coord, anno, label2thick = FALSE){
  # Wrapping data from annotation
  names(anno) <- mcols(anno)$name
  anno <- anno[seqnames(tx_coord)]
  exon_len <- mcols(anno)$blockSizes
  tx_len <- sum(exon_len)
  exon_start <- mcols(anno)$blockStarts + start(anno) - 1L
  exon_end <- exon_start + exon_len - 1L
  strand <- as.character(strand(anno))
  is_m <- strand == "-"

  # Wrapping data from coordination in transcript
  tx_coord$name <- as.character(seqnames(tx_coord))
  site_1 <- start(tx_coord)
  site_2 <- end(tx_coord)
  site_1[is_m] <- (tx_len - end(tx_coord) + 1)[is_m] # Flip the start site in transcript to plus strand
  site_2[is_m] <- (tx_len - start(tx_coord) + 1)[is_m] # Flip the end site in transcript to plus strand

  # Cumulative length of exon
  exon_cs <- cumsum(exon_len)

  # Finding where the stie falling at
  site_1_at <- unname(IntegerList(first_element(IRanges::which(site_1 <= exon_cs))))
  site_2_at <- unname(IntegerList(first_element(IRanges::which(site_2 <= exon_cs))))

  # Computing the genomic location
  genome_start <- exon_start[site_1_at] + (site_1 - (exon_cs[site_1_at] - exon_len[site_1_at])) - 1
  genome_start <- unlist(genome_start)
  genome_end <- exon_start[site_2_at] + (site_2 - (exon_cs[site_2_at] - exon_len[site_2_at])) - 1
  genome_end <- unlist(genome_end)
  out <- GRanges(seqnames(anno), IRanges(genome_start, genome_end),
          strand = strand)
  out$transcript_id <- anno$name
  return(out)
}

label_bed_thick <- function(bed, features, identifier = "transcript_id"){
  selected <- mcols(features)[[identifier]]
  if(any(duplicated(selected))){
    stop()
  }
  mcols(bed)[selected,]$thickStart <- start(features)
  mcols(bed)[selected,]$thickEnd <- end(features)
  return(bed)
}
