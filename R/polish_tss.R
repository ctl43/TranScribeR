#' @export
polish_tss <- function(ref_gtf, bed, cage, seq_info,
                       selected_fields = c("type","gene_id", "gene_name", "source", "transcript_id", "transcript_name", "exon_number", "exon_id")){
  ref_gtf <- import(ref_gtf)

  if(any(grepl("_", ref_gtf$transcript_id))){
    stop("Transcript_id should not contain '_', please remove them.")
  }

  ref_gtf <- ref_gtf[S4Vectors::match(seqnames(ref_gtf), seqnames(seq_info), nomatch = 0) > 0]
  seqlevels(ref_gtf) <- seqlevelsInUse(ref_gtf)
  seqinfo(ref_gtf) <- seq_info[seqlevels(ref_gtf)]

  tr <- ref_gtf[ref_gtf$type == "transcript"]
  names(tr) <- tr$transcript_id
  ex <- ref_gtf[ref_gtf$type == "exon"]
  ex <- get_first_exon(ex, label_only = TRUE)

  full_tss <- resize(tr, width = 1)
  tss <- unname(unique(full_tss))
  mcols(tss) <- NULL
  mcols(tss)$gene_id <- unique(full_tss)$gene_id
  tss$members <- CharacterList(split(full_tss$transcript_id, as.character(full_tss)))[as.character(tss)]

  # Importing cage seq data
  cage <- import_bed(cage, seqinfo = seq_info, select = "all", remove_thick = FALSE, use_name = FALSE)
  start(cage) <- cage$thickStart + 1
  end(cage) <- cage$thickStart + 1

  # Selecting confident TSS from predicted TSS
  legit <- unique(findOverlaps(cage, tss)@to)
  potential_ol <- findOverlaps(cage + 100, tss)
  potential_ol <- potential_ol[!potential_ol@to %in% legit]
  pr <- Pairs(cage[potential_ol@from], tss[potential_ol@to])
  tss_diff <- abs(start(pr@first) - start(pr@second))
  grp <- as.character(pr@second)
  grp <- factor(grp, levels = unique(grp))
  tss_diff <- IntegerList(split(tss_diff, grp))
  to_replace <- split(pr@first, grp)
  to_replace <- to_replace[tss_diff == min(tss_diff)]
  tmp_names <- names(to_replace)
  to_replace <- first_element(to_replace)
  names(to_replace) <- tmp_names
  to_replace <- unlist(to_replace)
  be_replaced <- unlist(first_element(split(pr@second, grp)))
  pr <- Pairs(to_replace, be_replaced)

  # Getting coverage
  reads <- import_bed(bed, select = "none", use_name = FALSE)
  read_start <- resize(reads, width = 1)
  read_start <- read_start[S4Vectors::match(seqnames(read_start), seqnames(seq_info), nomatch = 0) > 0]
  seqlevels(read_start) <- seqlevelsInUse(read_start)
  seqinfo(read_start) <- seq_info[seqlevels(read_start)]

  cage_cvg <- unlist(get_cvg(pr@first, read_start, append = FALSE, ignore_strand = FALSE))
  tss_cvg <- unlist(get_cvg(pr@second, read_start, append = FALSE, ignore_strand = FALSE))
  good_pr <- pr[(cage_cvg > 50) & (cage_cvg/tss_cvg) > 10]

  # Replacing the first exon that to be corrected
  new_tss <- rep(start(good_pr@first), lengths(good_pr@second$members))
  fe <- ex[ex$first_exon]
  mcols(fe)$first_exon <- NULL
  names(fe) <- fe$transcript_id
  be_replaced_exon <- fe[unlist(good_pr@second$members)]
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
  mcols(updated) <- mcols(updated)[names(mcols(updated))%in%selected_fields]
  updated <- updated[order(updated$transcript_id)]
  updated <- updated[order(updated$type)]
  updated <- sort(updated)
  updated$source <- "TranScribeR"
  return(updated)
}
