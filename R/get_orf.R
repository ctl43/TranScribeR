#' @export
#' @importFrom BiocGenerics paste
#' @importFrom Biostrings nchar
get_main_orf <- function(predicted_orf,
                         pfam_out = NULL,
                         blastp_out = NULL,
                         ref_start_codon = NULL,
                         anno = NULL,
                         min_aa = 50,
                         label_main_only = FALSE,
                         predict_nmd = TRUE){
  # To make sure when scoring metrics are same, the longer one is chosen
  orf <- predicted_orf
  orf <- orf[order(width(predicted$tx_coord), decreasing = TRUE), ]
  orf <- orf[order(seqnames(orf$tx_coord)), ]
  ref_start_ok <- FALSE
  if(!is.null(anno)){
    anno <- import_bed12(anno)
    orf$genome_coord <- tx2genome(orf$tx_coord, anno = anno)
    if(!is.null(ref_start_codon)){
      orf_start <- resize(orf$genome_coord, width = 1)
      orf$ref_start_supported <- ref_start_ok <- as.character(orf_start) %in% as.character(ref_start_codon)
    }
  }

  hmm_ok <- FALSE
  blastp_ok <- FALSE
  if(!is.null(blastp_out) & !is.null(pfam_out)){
    hmm <- fread(pfam_out, skip = 3, header = FALSE, fill = TRUE)
    hmm <- hmm[seq_len(nrow(hmm) - 10)]
    hmm[[7]] <- as.numeric(hmm[[7]])
    hmm <- hmm[hmm[[7]] < 1E-5]

    ## Getting result from blastp result
    blastp <- fread(blastp_out, header = FALSE)

    # Tidying up data
    hmm_ok <- names(orf$amino_acid) %in% hmm[[1]]
    blastp_ok <- names(orf$amino_acid) %in% blastp[[1]]
    orf_aa_id <- unique(names(orf$amino_acid))
    blastp_match <- CharacterList(split(blastp[[2]], factor(blastp[[1]], levels = orf_aa_id)))
    hmm_match <- CharacterList(split(hmm[[4]], factor(hmm[[1]], levels = orf_aa_id)))
    orf$pfam <- hmm_match[names(orf$amino_acid)]
    orf$blastp <- blastp_match[names(orf$amino_acid)]
  }

  # Getting longest predicted ORF within each transcript
  dt <- data.table(idx = seq_len(nrow(orf)), len = nchar(orf$amino_acid),
                   tx_ids = as.character(seqnames(orf$tx_coord)))
  selected <- dt[,.SD[which.max(len), ], by = tx_ids]
  longest_ok <- dt$idx %in% selected$idx
  orf$is_longest <- longest_ok

  # Computing metrics for selecting the best transcript d
  score <- longest_ok * 10 + ref_start_ok * 10 + hmm_ok + blastp_ok
  score[nchar(orf$amino_acid) < min_aa] <- 0 # ORF with length shorter than min length will have 0 score
  dt <- data.table(idx = seq_len(nrow(orf)), score = score,
                   tx_ids = as.character(seqnames(orf$tx_coord)))
  selected <- dt[,.SD[which.max(score), ], by = tx_ids]
  is_main_orf <- dt$idx %in% selected$idx
  is_main_orf[nchar(orf$amino_acid) < min_aa] <- FALSE
  is_main_orf[score <= 10] <- FALSE # longest + ribo supported | longest/ribo + hmm/blastp # Need to revise
  orf$main_orf <- is_main_orf

  if(predict_nmd){
    if(is.null(anno)){
      stop("NMD prediction requires transcript annotation files")
    }
    nmd_predicted <- predict_nmd(tx_orf = orf$tx_coord[orf$main_orf],
                                 anno = anno)
    orf$nmd <- NA
    orf$nmd[orf$main_orf] <- nmd_predicted
  }

  if(label_main_only){
    return(orf)
  }else{
    return(orf[is_main_orf, ])
  }
}
