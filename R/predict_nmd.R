#' @export
#' @importFrom rtracklayer import
#' @importFrom S4Vectors mcols 'mcols<-'
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics strand
predict_nmd <- function(genomic_cds = NULL, tx_orf = NULL, anno = NULL,
                        dist_2_lastjunc = 50, unique_id = c("transcript_id", "gene_name")){
  # Computing the length that a non-NMD transcript should be.
  # Transcript length - last exon length - dist_2_lastjunc (50)

  # Checking the annotation file
  if(is.character(anno)){
    is_file <- file.exists(anno)
    if(!is_file){
      stop("Please provide a valid file path.")
    }
  }else{
    is_file <- FALSE
  }

  if(is_file){
    anno_file <- anno
    file_type <- gsub(".*\\.","",basename(anno))
    if(file_type == "gtf"){
      anno <- gtf2bed12(anno)
    }else{
      if(file_type == "bed"){
        anno <- import_bed12(anno)
      }else{
        stop("An annotation file, in either bed12 or gtf format, should be provided.")
      }
    }
  }else{
    if(class(anno)=="GRanges"){
      is_bed12 <- all(c("blockSizes", "blockStarts") %in% names(mcols(anno)))
      is_gtf <- all(c("start_codon", "stop_codon", "exon") %in% (mcols(anno)$type))
      if((!is_bed12) & (!is_gtf)){
        stop("An annotation file, in either bed12 or gtf format, should be provided.")
      }
      if(is_gtf){
        anno <- gtf2bed12(anno, unique_id = unique_id)
      }
    }
  }

  last_bs <- rep(0, length(anno))
  is_p <- as.logical(strand(anno) == "+")
  last_bs[is_p] <- unlist(last_element(anno[is_p]$blockSizes), use.names = FALSE)
  last_bs[!is_p] <- unlist(first_element(anno[!is_p]$blockSizes), use.names = FALSE)
  transcript_len <- sum(anno$blockSizes)
  no_last_len <- transcript_len - last_bs
  names(no_last_len) <- names(anno)

  if(!is.null(tx_orf)){
    # Must be CDS
    coord <- tx_orf
    names(coord) <- seqnames(coord)
    is_nmd <- rep(NA, length(coord))
    names(is_nmd) <- names(coord)
  }

  if(!is.null(genomic_cds) & is.null(tx_orf)){
    if(class(genomic_cds)=="character"){
      if(genomic_cds == anno_file){
        genomic_cds <- anno
      }
    }
    if(class(genomic_cds)=="GRanges"){
      is_bed12 <- all(c("blockSizes", "blockStarts", "thickStart", "thickEnd") %in% names(mcols(genomic_cds)))
      is_gtf <- all(c("start_codon", "stop_codon", "exon") %in% (mcols(genomic_cds)$type))
      if(is_gtf){
        genomic_cds <- gtf2bed12(genomic_cds, unique_id = unique_id)
        coord <- genome2tx(anno = genomic_cds)
      }
      if(is_bed12){
        coord <- genome2tx(anno = genomic_cds)
      }
      if(!is_gtf & !is_bed12){
        coord <- genome2tx(target = genomic_cds, anno = anno)
      }

    }
    # Can contain annotation that may not code for protein
    names(coord) <- seqnames(coord)
    is_nmd <- rep(NA, length(genomic_cds))
    names(is_nmd) <- names(genomic_cds)
    coord <- coord[start(coord)!=0]
  }

  if(any(duplicated(as.character(seqnames(coord))))){
    stop("Each transcript can only have one main ORF.")
  }

  selected_is_nmd <- (end(coord) - (no_last_len[names(coord)] - dist_2_lastjunc)) < 0
  is_nmd[names(coord)] <- selected_is_nmd
  return(is_nmd)
}
