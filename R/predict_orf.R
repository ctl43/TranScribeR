#'@export
find_codon <- function(codon, seq){
  putative <- gregexpr(codon, as.character(seq))
  putative <- IntegerList(putative)
  putative[putative < 0] <- 0
  grp <- get_list_grp(putative)
  seqnames <- rep(names(seq), lengths(putative))
  putative <- GRanges(factor(seqnames, unique(names(seq))), IRanges(unlist(putative), width = 1))
  putative <- putative[start(putative) != 0]
  putative$codon <- codon
  putative$phase <- start(putative) %% 3
  return(putative)
}

#'@export
predict_orf <- function(fa, seq = NULL, start_codon = "ATG", starts = NULL,
                        stop_codon = c("TAG", "TAA", "TGA"),
                        min_aa = 30, ref_prot_fa = NULL,
                        assign_prot_id = TRUE, with_stop = FALSE){
  # ORF predicted in here must have start and stop codon.
  if(is.null(seq)){
    seq <- read_fasta(fa)
  }
  seq <- as.character(seq)

  if(!is.null(start_codon)){
    starts <- lapply(start_codon, find_codon, seq = seq)
    starts <- unlist(GRangesList(starts), use.names = FALSE)
  }
  seq <- seq[seqlevels(starts)]
  phased_start <- split(starts, factor(starts$phase, levels = 0L:2L))
  phased_start <- lapply(phased_start, function(x)split(start(x), seqnames(x)))
  phased_start <- lapply(phased_start, function(x)x - 1L)
  phased_predicted <- lapply(phased_start, cxx_find_orf, seqs = seq, stop_codon = stop_codon)
  phased_predicted <- DataFrame(lapply(phased_predicted, function(x)IntegerList(x) + 1L))
  colnames(phased_predicted) <- c("frame0", "frame1", "frame2")

  phased_predicted <- lapply(phased_predicted, function(x){
    seqnames <- rep(names(seq), lengths(x)/2)
    y <- matrix(unlist(x, use.names = FALSE), ,2, byrow = TRUE)
    GRanges(seqnames, IRanges(y[,1], y[,2]))
  })
  predicted <- GRangesList(phased_predicted)
  phase <- Rle(rep(0:2, lengths(predicted)))
  predicted <- unlist(predicted, use.names = FALSE)
  extracted <- DNAStringSet(subseq(seq[as.character(seqnames(predicted))], start(predicted), end(predicted)))
  aa <- translate(extracted, no.init.codon = TRUE)
  aa <- sub("\\*$", "", aa)
  aa <- AAStringSet(aa)
  first_codon <- as.character(subseq(extracted, start = 1, end = 3))
  last_codon <- unname(as.character(subseq(extracted, start = lengths(extracted) - 2, end = lengths(extracted))))
  ids <- do.call(rbind, strsplit(as.character(seqnames(predicted)), "_"))
  out <- DataFrame(transcript_id = ids[,1], gene_id = ids[,2],
                   input_length = nchar(seq[as.character(seqnames(predicted))]),
                   tx_coord = predicted, phase = phase, cDNA = extracted, amino_acid = aa,
                   first_codon = first_codon, last_codon = last_codon, has_stop = match(last_codon, stop_codon, nomatch = 0) > 0)

  if(assign_prot_id){
    # Assigning protein id
    out <- out[order(names(out$amino_acid)), ]
    unique_prot <- out[!duplicated(out$amino_acid), ]
    protein_id <- as.character(unique_prot$tx_coord)
    names(protein_id) <- as.character(unique_prot$amino_acid)

    if(!is.null(ref_prot_fa)){
      # Referencing to the annotated protein id
      ref_pep <- read_fasta(ref_prot_fa)
      ref_pep <- ref_pep[!duplicated(ref_pep)]
      ref_protein_id <- names(ref_pep)
      names(ref_protein_id) <- as.character(ref_pep)
      cm_pep <- intersect(names(protein_id), as.character(ref_pep))
      protein_id[cm_pep] <- ref_protein_id[cm_pep]
    }
    names(out$amino_acid) <- unname(protein_id[as.character(out$amino_acid)])
  }

  if(with_stop){
    out <- out[out$has_stop, ]
  }
  if(min_aa > 0){
    out <- out[lengths(out$amino_acid) >= min_aa,]
  }
  return(out)
}
