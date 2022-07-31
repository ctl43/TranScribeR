#' @export
#' @importFrom BiocGenerics paste
read_fasta <- function(x){
  # y <- read.fasta(x)
  y <- fread(x, sep = "\t", header = FALSE)[[1]]
  is_id <- grepl("^>", y)
  grp <- cumsum(is_id)
  seq_id <- y[is_id]
  seqs <- y[!is_id]
  seq_grp <- grp[!is_id]
  seqs <- paste(CharacterList(split(seqs, factor(seq_grp, levels = seq_len(sum(is_id))))), collapse = "")
  seqs <- unlist(seqs)
  # anno <- sapply(y, function(x)attr(x, "Annot"))
  seqs <- BStringSet(toupper(seqs))
  names(seqs)  <- sub(">","",seq_id)
  # mcols(seqs)$annotation <- anno
  seqs
}

#' @export
write_fasta <- function(x, out){
  if(class(x)=="BStringSet"){
    y <- as.character(x)
    names(y) <- names(x)
    x <- y
  }
  fa_names <- paste0(">", names(x))
  write(as.character(rbind(fa_names, x)), out)
}

#' @export
write_fastq <- function(id, seq, qual, out){
  result <- as.character(rbind(paste0("@",id), seq, paste0("+", id), qual))
  write(result, out)
}

#' @export
#' @importFrom BiocGenerics start
get_frame <- function(cds){
  (start(cds) - 1)%%3
}

#' @export
first_element <- function(x, invert = FALSE){
  if(length(x) == 0){
    return(x)
  }
  grp <- factor(rep(seq_along(x), lengths(x)), levels = seq_along(x))
  flat <- unlist(x, use.names = FALSE)
  is_dup <- duplicated(grp)
  if(invert){
    return(split(flat[is_dup], grp[is_dup]))
  }else{
    return(split(flat[!is_dup], grp[!is_dup]))
  }
}

#' @export
last_element <- function(x, invert = FALSE){
  if(length(x) == 0){
    return(x)
  }
  grp <- factor(rep(seq_along(x), lengths(x)), levels = seq_along(x))
  flat <- unlist(x, use.names = FALSE)
  rev_grp <- rev(grp)
  rev_flat <- rev(flat)
  is_dup <- duplicated(rev_grp)
  if(invert){
    return(split(rev(rev_flat[is_dup]), rev(rev_grp[is_dup])))
  }else{
    return(split(rev(rev_flat[!is_dup]), rev(rev_grp[!is_dup])))
  }
}

#' @export
get_list_grp <- function(x, as_factor = TRUE){
  tot <- seq_along(x)
  if(as_factor){
    tot <- factor(tot)
  }
  return(rep(tot, lengths(x)))
}

#' @export
#' @importFrom S4Vectors split
merge_granges <- function(x, tol, ignore_strand = TRUE){
  is_grlist <- class(x)=="CompressedGRangesList"
  if(is_grlist){
    grp <- get_list_grp(x, as_factor = TRUE)
    gr <- unlist(x)
  }else{
    gr <- x
    grp <- rep(1, length(gr), as_factor = TRUE)
  }
  chrom <- as.integer(seqnames(gr))
  start <- as.integer(start(gr))
  end <- as.integer(end(gr))

  if(ignore_strand){
    strand <- rep(1, length(gr), as_factor = TRUE)
    gr <- unstrand(gr)
  }else{
    strand <- as.integer(factor(strand(gr), levels = c("+", "-", "*")))
  }

  o <- order(grp, strand, chrom, start, end) # MUST BE SORTED LIKE THIS
  idx <- rep(0, length(gr))
  idx[o] <- cxx_merge_ranges(chrom[o], start[o], end[o], grp = grp[o], strand = strand[o], tol = tol)
  y <- IntegerList(split(idx, grp))
  members <- split(gr, unlist(idx))
  megred <- unlist(range(members))
  out_grp <- y - unname(cumsum(c(0, head(lengths(unique(y)), -1)))) # converting the idx
  out_grp <- IntegerList(out_grp)

  if(is_grlist){
    merged_grp <- unique(y)
    out_merged <- split(megred, get_list_grp(merged_grp, as_factor = TRUE))
    names(out_grp) <- names(out_merged) <- names(x)
    return(list(idx = out_grp, regions = out_merged))
  }else{
    return(list(idx = unlist(out_grp, use.names = FALSE), regions = megred))
  }
}

#' @export
#' @importFrom S4Vectors DataFrame mcols "mcols<-"
#' @importFrom data.table fread
#' @importFrom IRanges IntegerList IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqinfo
import_bed12 <- function(x,
                       remove_thick = TRUE,
                       select = "all",
                       use_name = TRUE,
                       seqinfo = NULL){
  touch <- fread(x, nrows = 1, sep = "\t")
  n <- ncol(touch)
  col_names <- c("chrom","start","end",
                 "name", "score","strand",
                 "thickStart", "thickEnd", "itemRgb",
                 "blockCount", "blockSizes", "blockStarts")
  fields_class <- c("character", "integer", "integer",
                    "character", "character", "character",
                    "integer", "integer", "character",
                    "character", "character", "character")
  bed <- DataFrame(fread(x, colClasses = fields_class, header = FALSE, sep = "\t"))
  colnames(bed) <- col_names

  bed[["blockSizes"]] <- IntegerList(strsplit(bed[["blockSizes"]], ","))
  bed[["blockStarts"]] <- IntegerList(strsplit(bed[["blockStarts"]], ","))
  bed$thickStart <- bed$thickStart + 1
  bed$blockStarts <- bed$blockStarts + 1
  bed$start <- bed$start + 1
  if(remove_thick){
    bed$thickStart <- 0
    bed$thickEnd <- 0
  }
  gr <- GRanges(bed$chrom, IRanges(bed$start, bed$end), strand = bed$strand, seqinfo = seqinfo)
  bed$name <- factor(bed$name, levels = unique(bed$name))
  if(use_name){
    names(gr) <- bed$name
  }
  mcols(gr) <- bed[-c(1,2,3,6)]
  return(gr)
}

#' @export
export_bed12 <- function(bed, out){
  df <- as.data.frame(bed)
  df <- DataFrame(df)
  colnames(df)[colnames(df)=="seqnames"] <- "chrom"
  # info <- mcols(bed)
  bed12field <- c("chrom","start","end", "name", "score","strand",
                  "thickStart", "thickEnd", "itemRgb", "blockCount",
                  "blockSizes", "blockStarts")
  if(!all(bed12field%in%names(df))){
    stop()
  }else{
    df <- df[, bed12field]
  }
  df$start <- df$start - 1L
  df$thickStart <- df$thickStart - 1L
  df$blockSizes <- paste0(paste(IntegerList(df$blockSizes), collapse = ","), ",")
  df$blockStarts <- paste0(paste(IntegerList(df$blockStarts) - 1L, collapse = ","), ",")
  write.table(df, out, quote=FALSE,
              col.names=FALSE, row.names = FALSE,
              sep = "\t")

}


#' @export
export_unique_pep <- function(orf, out, identifier = "protein_id"){
  if(!(identifier %in% colnames(orf))){
    stop("The selected data are not found.")
  }
  coding <- orf[!is.na(orf$protein_id), ]
  coding <- coding[!duplicated(coding$protein_id),]
  protein_id <- coding$protein_id
  aa_seq <- coding$aa_seq
  names(aa_seq) <- protein_id
  write_fasta(aa_seq, out)
}

#'@export
unlist_s4list <- function(y){
  # Work but clumsy (improve this later)
  what <- class(y[[1]])
  tmp <- do.call(c, lapply(y, function(x){
    split(unlist(x, use.names = FALSE), get_list_grp(x))
  }))
  unname(get(sub("Compressed","",what))(tmp))
}

#' @export
#' @importFrom S4Vectors split
unstranded_get_cvg <- function(gr, cvg, append = TRUE, ignore_strand = TRUE){
  if(!all(seqlevels(gr) %in% names(cvg))){
    stop()
  }
  cvg <- cvg[seqlevels(gr)]
  o <- order(unlist(split(seq_along(gr), seqnames(gr)), use.names = FALSE))
  gr$idx <- seq_along(gr)
  gr <- split(gr, seqnames(gr))
  region_cvg <- RleViewsList(mapply(function(x, y)Views(y, start(x), end(x)), x = gr, y = cvg))
  region_cvg <- lapply(region_cvg, IntegerList)
  region_cvg <- unlist_s4list(region_cvg)
  if(append){
    gr <- unlist(gr, use.names = FALSE)
    gr$cvg <- region_cvg
    return(gr[o])
  }else{
    return(region_cvg[o])
  }
}

#' @export
#' @importFrom S4Vectors split
get_cvg <- function(gr, cvg_gr, append = FALSE, ignore_strand = FALSE){
  if(!ignore_strand){
    cvg_gr <- split(cvg_gr, strand(cvg_gr))
    cvg <- lapply(cvg_gr, coverage)
    s_o <- order(unlist(split(seq_along(gr), strand(gr)), use.names = FALSE))
    gr <- split(gr, strand(gr))
    cvg_out <- mapply(unstranded_get_cvg, gr = gr, cvg = cvg[names(gr)], append = FALSE, SIMPLIFY = FALSE)
    cvg_out <- unlist_s4list(cvg_out)
    gr <- unlist(gr)
    if(append){
      gr$cvg <- cvg_out
      return(gr[s_o])
    }else{
      return(cvg_out[s_o])
    }
  }else{
    cvg <- coverage(cvg_gr)
    return(unstranded_get_cvg(gr, cvg, append = append))
  }
}

#'@export
#'@importFrom rtracklayer import
gtf2bed12 <- function(gtf, unique_id = c("transcript_id", "gene_name"), label_cds = TRUE){
  # when label cds is on, transcripts with multiple or with either start/stop codon will not be labeled.
  if(class(gtf)!="GRanges"){
    gtf <- import(gtf)
  }
  gtf <- unname(gtf)
  exon <- gtf[gtf$type=="exon"]
  exon <- sort(exon)
  exon <- split(exon, Reduce(function(...)paste(..., sep = "_"), as.list(mcols(exon[, unique_id]))))
  exon <- exon[order(names(exon))]
  gr <- unlist(range(exon))
  gr$name <- names(gr)
  gr$thickEnd <- gr$thickStart <- gr$score <- 0L
  gr$itemRgb <- "0"
  gr$blockCount <- lengths(exon)
  gr$blockSizes <- width(exon)
  gr$blockStarts <- unname(start(exon)) - unname(min(start(exon))) + 1

  # If there are start and stop codon
  start_codon <- gtf[gtf$type=="start_codon"]
  start_codon <- start_codon[width(start_codon) == 3]
  stop_codon <- gtf[gtf$type=="stop_codon"]
  stop_codon <- stop_codon[width(stop_codon) == 3]

  if(length(start_codon) > 0 & length(stop_codon) > 0 & label_cds){
    start_codon <- split(start_codon, Reduce(function(...)paste(..., sep = "_"), as.list(mcols(start_codon[, unique_id]))))
    stop_codon <- split(stop_codon, Reduce(function(...)paste(..., sep = "_"), as.list(mcols(stop_codon[, unique_id]))))
    if(any(lengths(start_codon) > 1) | any(lengths(stop_codon) > 1)){
      stop("Each gene should not have more than one start codon/stop codon")
    }else{
      start_codon <- unlist(start_codon)
      start_codon <- resize(start_codon, width = 1)
      start_codon <- start_codon[order(names(start_codon))]
      stop_codon <- unlist(stop_codon)
      stop_codon <- resize(stop_codon, width = 1, fix = "end")
      stop_codon <- stop_codon[order(names(stop_codon))]
    }

    cm <- intersect(names(start_codon), names(stop_codon))
    start_codon <- start_codon[cm]
    stop_codon <- stop_codon[cm]
    is_p <- names(gr) %in% cm & strand(gr) == "+"
    is_m <- names(gr) %in% cm & strand(gr) == "-"
    gr[is_p]$thickStart <- start(start_codon[strand(start_codon) == "+"])
    gr[is_p]$thickEnd <- start(stop_codon[strand(stop_codon) == "+"])
    gr[is_m]$thickEnd <- start(start_codon[strand(start_codon) == "-"])
    gr[is_m]$thickStart <- start(stop_codon[strand(stop_codon) == "-"])
  }
  return(gr)
}

#'@export
bed122gtf <- function(bed, name_sep = "_"){
  # Add an assume_cds = FALSE option in the future
  if(class(bed)=="character"){
    bed <- import_bed12(bed)
  }else{
    if(class(bed)!="GRanges"){
      stop()
    }
  }
  bed$blockCount <- as.integer(bed$blockCount)
  exon_chr <- rep(seqnames(bed), bed$blockCount)
  exon_strand <- rep(strand(bed), bed$blockCount)
  gene <- paste(CharacterList(last_element(strsplit(as.character(bed$name), name_sep))), collapse = "_")
  gene <- unname(gene)
  exon_gene <- rep(gene, bed$blockCount)
  tx <- paste(CharacterList(last_element(strsplit(as.character(bed$name), name_sep), invert = TRUE)), collapse = "_")
  tx <- unname(tx)
  exon_tx <- rep(tx, bed$blockCount)
  exon_start <- unlist(bed$blockStarts) + rep(start(bed), bed$blockCount)
  exon_end <- exon_start + unlist(bed$blockSizes) - 1
  gr <- GRanges(exon_chr, IRanges(unlist(exon_start), unlist(exon_end)), strand = exon_strand)
  mcols(bed) <- NULL
  bed$source <- gr$source <- "TranScribeR"
  bed$type <- "transcript"
  gr$type <- "exon"
  bed$score <- gr$score <- NA
  bed$phase <- gr$phase <- as.integer(NA)
  gr$gene_id <- exon_gene
  bed$gene_id <- gene
  gr$transcript_id <- exon_tx
  bed$transcript_id <- tx
  bed$exon_number <- gr$exon_number <- as.integer(NA)
  gr <- c(bed, gr)
  gr <- unname(gr)
  gr <- sort(gr)
  return(gr)
}

#'@export
characterise <- function(x, list_only = FALSE){
  x <- DataFrame(x)
  classes <- sapply(x, class)
  x[,classes=="list"] <- DataFrame(lapply(x[,classes=="list", drop = FALSE], CharacterList, collapse = ","))
  classes <- sapply(x, class)
  is_s4list <- grepl("List", classes)
  x[,is_s4list] <- DataFrame(lapply(x[,is_s4list, drop = FALSE], paste, collapse = ","))
  if(!list_only){
    x[,!is_s4list] <- DataFrame(lapply(x[,!is_s4list, drop = FALSE], as.character))
  }
  x
}

#'@export
annotate_morf <- function(x, anno, out){
  bed <- import_bed12(anno, use_name = TRUE, remove_thick = TRUE)
  morf <- orf[orf$main_orf, ]
  bed[as.character(seqnames(morf$tx_coord))]$thickStart <- start(morf$genome_coord)
  bed[as.character(seqnames(morf$tx_coord))]$thickEnd <- end(morf$genome_coord)
  export_bed12(bed, out)
}

#'@export
append_noncoding <- function(morf, fa){
  # fa <- "flair_out/collapsed/using_all_reads/assigned_collapsed.isoforms.fa"
  seq <- read_fasta(fa)
  seq <- seq[!names(seq) %in% as.character(seqnames(morf$tx_coord))]
  non_coding_set <- names(seq)
  ids <- do.call(rbind, strsplit(non_coding_set, "_"))
  colnames(ids) <- c("transcript_id", "gene_id")
  input_length <- nchar(seq[non_coding_set])

  out <- DataFrame(DataFrame(ids),
                   input_length = input_length,
                   tx_coord = GRanges(non_coding_set, IRanges(0)))

  for(z in colnames(morf)[5:ncol(morf)]){
    i <- class(morf[[z]])
    if(i == "character"){
      out[[z]] <- ""
    }
    if(i == "logical"){
      out[[z]] <- FALSE
    }
    if(i == "CompressedGRangesList"){
      out[[z]] <- GRangesList(GRanges())
    }
    if(i == "Rle"){
      out[[z]] <- Rle(NA)
    }
    if(i == "DNAStringSet"){
      out[[z]] <- DNAStringSet("N")
    }
    if(i == "AAStringSet"){
      out[[z]] <- AAStringSet("N")
    }
    if(i == "CompressedCharacterList"){
      out[[z]] <- CharacterList(character())
    }
    if(i == "GRanges"){
      out[[z]] <- GRanges(" ", IRanges(0))
    }
  }

  out <- DataFrame(rbind(morf[,1:3], out[,1:3]),
                   tx_coord = suppressWarnings(c(morf$tx_coord, out$tx_coord)),
                   suppressWarnings(rbind(morf[,5:ncol(morf)], out[, 5:ncol(out)])),
                   coding = c(rep(TRUE, nrow(morf)), rep(FALSE, nrow(out))))

  return(out)
}

#'@export
append_gene <- function(anno){
  tx <- anno[anno$type=="transcript"]
  gene <- range(split(tx, tx$gene_id))
  if(any(lengths(gene)!=1)){
    stop("Some genes are not on the same chromosome or strand.")
  }
  gene_id <- names(gene)
  gene <- unlist(gene)
  gene$source <- "TranScribeR"
  gene$type <- "gene"
  gene$gene_id <- gene_id
  c(anno, unname(gene))
}

#'@export
append_CDS <- function(anno, tx_coord, phase){
  anno_bed <- gtf2bed12(anno, unique_id = c("transcript_id", "gene_id"))
  genome_coord <- tx2genome(tx_coord, anno = anno_bed)
  grrl <- split(genome_coord, genome_coord$transcript_id)
  exon <- anno[anno$type == "exon"]
  exon <- sort(exon)
  exon <- split(exon, paste0(exon$transcript_id, "_", exon$gene_id))
  exon <- exon[names(grrl)]
  cds <- restrict_list(exon, grrl)
  ids <- names(exon)
  ids <- rep(ids, lengths(cds))
  phase <- rep(phase, lengths(cds))
  cds <- unlist(cds, use.names = FALSE)
  cds$type <- "CDS"
  ids <- do.call(rbind, strsplit(ids, "_"))
  cds$gene_id <- ids[,2]
  cds$transcript_id <- ids[,1]
  cds$source <- "TranScribeR"
  cds$phase <- phase
  return(c(anno, cds))
}

#'@export
restrict_list <- function(grl, grrl){
  # Need to be polished
  if(any(lengths(grrl) != 1)){
    stop()
  }
  grl_grp <- get_list_grp(grl)
  grl <- unlist(grl, use.names = FALSE)
  mcols(grl) <- NULL
  grl <- split(grl, grl_grp)

  # Filter out uncommon seqnames
  grl <- grl[seqnames(grl)==seqnames(grrl)]

  # Filter out not in range
  grl <- grl[!(end(grl) < start(grrl) | start(grl) > end(grrl))]

  # Getting those not out of range and save it
  ok <- (!start(grl) < start(grrl) & !end(grl) > end(grrl))
  ok_gr <- grl[ok]
  ok_grp <- get_list_grp(ok_gr)
  ok_gr <- unlist(ok_gr, use.names = FALSE)

  # Getting those out of range
  ## Left first
  grl <- grl[!ok]
  too_left <- start(grl) < start(grrl)
  too_left_gr <- grl[too_left]
  too_left_grp <- get_list_grp(too_left_gr)

  # retrict them into correct range
  left_sub <- rep(start(grrl), lengths(too_left_gr))
  left_sub <- unlist(left_sub, use.names = FALSE)
  too_left_gr <- unlist(too_left_gr, use.names = FALSE)
  start(too_left_gr) <- left_sub

  non_too_left_gr <- grl[!too_left]
  non_too_left_grp <- get_list_grp(non_too_left_gr)
  non_too_left_gr <- unlist(non_too_left_gr, use.names = FALSE)


  grl <- split(c(too_left_gr, non_too_left_gr),
               factor(c(too_left_grp, non_too_left_grp),
                      levels = levels(too_left_grp)))

  # Repeat the same procedure for right
  too_right <- end(grl) > end(grrl)
  too_right_gr <- grl[too_right]
  too_right_grp <- get_list_grp(too_right_gr)

  # retrict them into correct range
  right_sub <- rep(end(grrl), lengths(too_right_gr))
  right_sub <- unlist(right_sub, use.names = FALSE)
  too_right_gr <- unlist(too_right_gr, use.names = FALSE)
  end(too_right_gr) <- right_sub
  non_too_right_gr <- grl[!too_right]
  non_too_right_grp <- get_list_grp(non_too_right_gr)
  non_too_right_gr <- unlist(non_too_right_gr, use.names = FALSE)

  grl <- split(c(too_right_gr, non_too_right_gr, ok_gr),
               factor(c(too_right_grp, non_too_right_grp, ok_grp),
                      levels = levels(ok_grp)))

  final_grp <- get_list_grp(grl)
  flat <- unlist(grl, use.names = FALSE)
  new_order <- order(flat)
  flat <- flat[new_order]
  final_grp <- final_grp[new_order]
  split(flat, final_grp)
}
