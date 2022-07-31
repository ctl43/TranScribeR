#' @export
#' @importFrom S4Vectors match
#' @importFrom rtracklayer import
merge_tss <- function(ref_gtf, seq_info, bed, tol = 20, min_cvg = 10,
                      selected_fields = c("type","gene_id", "gene_name", "source", "transcript_id", "transcript_name", "exon_number", "exon_id")){

  # Importing annotation
  ref_gtf <- import(ref_gtf)
  if(any(grepl("_", ref_gtf$transcript_id))){
    stop("Transcript_id should not contain '_', please remove them.")
  }
  ref_gtf <- ref_gtf[S4Vectors::match(seqnames(ref_gtf), seqnames(seq_info), nomatch = 0) > 0]
  seqlevels(ref_gtf) <- seqlevelsInUse(ref_gtf)
  seqinfo(ref_gtf) <- seq_info[seqlevels(ref_gtf)]
  tr <- ref_gtf[ref_gtf$type == "transcript"]
  names(tr) <- tr$transcript_id
  tss <- resize(tr, width = 1)
  ex <- ref_gtf[ref_gtf$type == "exon"]
  ex <- get_first_exon(ex, label_only = TRUE)

  # Getting the coverage
  reads <- import_bed(bed, select = "none", use_name = FALSE)
  read_start <- resize(reads, width = 1)
  read_start <- read_start[S4Vectors::match(seqnames(read_start), seqnames(seq_info), nomatch = 0) > 0]
  seqlevels(read_start) <- seqlevelsInUse(read_start)
  seqinfo(read_start) <- seq_info[seqlevels(read_start)]
  tss$cvg <- unlist(get_cvg(tss, read_start, append = FALSE, ignore_strand = FALSE))

  # Merging TSS if they are close enough and choose the one with highest coverage
  tss <- split(tss, tss$gene_id)
  merged_tss <- merge_granges(tss, tol = tol, ignore_strand = FALSE)
  tmp_grp <- paste0(rep(names(merged_tss[[1]]), lengths(merged_tss[[1]])),
                    "_", unlist(merged_tss[[1]]))
  tss <- unlist(tss)
  tss_members <- split(tss$transcript_id, tmp_grp)
  cvg_dt <- data.table(chrom = as.character(seqnames(tss)),
                       start = as.integer(start(tss)),
                       strand = as.character(strand(tss)),
                       cvg = tss$cvg, grp = tmp_grp)
  selected_tss <- cvg_dt[, .SD[which.max(cvg),], by = grp]

  # For TSS without any supporting reads, the TSS will not be corrected
  selected_tss <- selected_tss[selected_tss$cvg >= min_cvg,]
  selected_members <- tss_members[selected_tss$grp]
  new_tss <- rep(selected_tss[[3]], lengths(selected_members))
  fe <- ex[ex$first_exon]
  names(fe) <- fe$transcript_id
  be_replaced_exon <- fe[unlist(selected_members)]
  is_p <- as.logical(strand(be_replaced_exon)=="+")

  # Updating the First exon annotation
  # For positive strand
  p <- be_replaced_exon[is_p]
  p_tss <- new_tss[is_p]
  p_discard <- p_tss > end(p)
  p_tss <- p_tss[!p_discard]
  p <- p[!p_discard]
  start(p) <- p_tss

  # For minus strand
  m <- be_replaced_exon[!is_p]
  m_tss <- new_tss[!is_p]
  m_discard <- m_tss < start(m)
  m_tss <- m_tss[!m_discard]
  m <- m[!m_discard]
  end(m) <- m_tss
  new_fe <- c(p, m)
  fe[names(new_fe)] <- new_fe

  # Updating transcript annotation
  start(tr[p$transcript_id]) <- start(p)
  end(tr[m$transcript_id]) <- end(m)

  # Reconstruct gene annotation
  ge <- range(split(tr, tr$gene_id))
  id2names <- tr[!duplicated(tr$gene_id)]$gene_name
  names(id2names) <- tr[!duplicated(tr$gene_id)]$gene_id
  gene_ids <- rep(names(ge), lengths(ge))
  ge <- unlist(ge)
  ge$gene_id <- gene_ids
  ge$gene_name <- id2names[gene_ids]
  ge$type <- "gene"

  # Combining the annotation
  updated <- unname(c(fe, ex[!ex$first_exon], tr, ge))
  mcols(updated) <- mcols(updated)[names(mcols(updated)) %in% selected_fields]
  updated <- updated[order(updated$transcript_id)]
  updated <- updated[order(updated$type)]
  updated <- sort(updated)
  updated$source <- "TranScribeR"
  mcols(updated)$first_exon <- NULL
  updated
}
