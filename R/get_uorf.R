#' @export
get_uorf <- function(tx_coord, phase, fa, anno, ref_start,
                     start_codons = c("ATG", "TTG","ACG", "CTG","GTG"),
                     stop_codons = c("TAG", "TAA", "TGA"), min_aa = 2){
  phase <- as.integer(phase)
  # tx_coord <- morf$tx_coord
  names(phase) <- seqnames(tx_coord)

  # Limiting start codons at five prime UTR
  fiveputr <- get_5putr(tx_coord, min_len = 9)
  fa <- read_fasta(fa)
  fivep_seq <- subseq(fa[seqnames(fiveputr)], fiveputr@ranges)
  starts <- lapply(start_codons, find_codon, seq = fivep_seq)
  starts <- unlist(GRangesList(starts), use.names = FALSE)
  starts <- starts[starts$phase != phase[as.character(seqnames(starts))]] # Not the same phase as main ORFs

  # Selected start codon supported by reference start site
  if(!is.null(ref_start)){
    anno <- import_bed12(anno)
    genomic <- tx2genome(starts, anno = anno)
    supp <- starts[S4Vectors::match(genomic, ref_start, nomatch = 0) > 0]
    supp <- predict_orf(seq = fa, starts = supp, start_codon = NULL,
                        stop_codon = stop_codons, min_aa = min_aa,
                        assign_prot_id = FALSE)
    supp$ref_supported <- TRUE
    unsupp <- starts[!S4Vectors::match(genomic, ref_start, nomatch = 0) > 0]
    unsupp <- predict_orf(seq = fa, starts = unsupp, start_codon = NULL,
                          stop_codon = stop_codons, min_aa = min_aa,
                          assign_prot_id = FALSE)
    unsupp$ref_supported <- FALSE
    supp <- S4Vectors::split(supp, factor(supp$phase, levels = 0:2))
    unsupp <- S4Vectors::split(unsupp, factor(unsupp$phase, levels = 0:2))

    # Selecting unsupported uORF only when there is no overlapping uORF supported by reference start site
    suppressWarnings(unsupp_only <- mapply(function(x, y){
      ol <- findOverlaps(x$tx_coord, y$tx_coord)
      y[!seq_len(nrow(y)) %in% ol@to, ]
    }, x = supp, y = unsupp, SIMPLIFY = FALSE))
    suppressWarnings(unsupp_only <- do.call(rbind, unsupp_only))
    supp <- unlist(supp, use.names = FALSE)
    uorf <- suppressWarnings(rbind(supp, unsupp_only))

  }else{
    # Basically it takes all ORF in 5' UTR
    uorf <- predict_orf(seq = fa, starts = starts, start_codon = NULL,
                stop_codon = stop_codons, min_aa = min_aa,
                assign_prot_id = FALSE)
    uorf$ref_supported <- FALSE

  }
  extra <- seqlevels(tx_coord)[!seqlevels(tx_coord)%in%seqlevels(uorf$tx_coord)]
  seqlevels(uorf$tx_coord) <- c(seqlevels(uorf$tx_coord), extra)
  return(uorf)
}
