#' @export
get_junctions <- function(x, gene_id = "gene_id", transcript_id = "transcript_id"){

  # flair <- flair[lengths(flair$blockSizes)!=1]
  # left_site <- start(flair) + flair$blockStarts - 1
  # right_site <- left_site + flair$blockSizes
  # left_j <- IntegerList(first_element(left_site, invert = TRUE))
  # right_j <- IntegerList(last_element(right_site, invert = TRUE))
  # j_loc <- IntegerList(split(c(unlist(left_j, use.names = FALSE),
  #                              unlist(right_j, use.names = FALSE)),
  #                            c(get_list_grp(left_j), get_list_grp(right_j))))

  exon <- x[x$type == "exon"]
  junctions <- sort(exon)
  junctions <- split((junctions), junctions$transcript_id)
  junctions <- junctions[lengths(junctions) != 1]

  chr <- seqnames(junctions)
  chr <- first_element(chr, invert = TRUE)
  chr <- unlist(chr)

  strand <- strand(junctions)
  strand <- first_element(strand, invert = TRUE)
  strand <- unlist(strand)

  j_start <- IntegerList(first_element(start(junctions), invert = TRUE)) + 1
  j_end <- IntegerList(last_element(end(junctions), invert = TRUE)) - 1
  # Faster than sort(IntegerList(mapply(c, j_start, j_end))) but look clumsier
  junctions_loc <- IntegerList(split(c(unlist(j_start, use.names = FALSE),
                                       unlist(j_end, use.names = FALSE)),
                                     c(get_list_grp(j_start), get_list_grp(j_end))))
  junctions_loc <- sort(junctions_loc)
  junctions_loc <- unlist(junctions_loc)
  junc_start <- junctions_loc[seq(1, length(junctions_loc), by = 2)]
  junc_end <- junctions_loc[seq(2, length(junctions_loc), by = 2)]
  junc <- GRanges(chr, IRanges(junc_start, junc_end), strand = strand)

  grp <- rep(seq_along(junctions), lengths(junctions))
  junc_grp <- unlist(first_element(split(grp, grp), invert = TRUE))
  junc <- split(junc, junc_grp)
  names(junc) <- names(junctions)

  junc_genes <- mcols(x)[[gene_id]][!duplicated(x$transcript_id)]
  names(junc_genes) <- x$transcript_id[!duplicated(x$transcript_id)]
  mcols(junc)$genes <-  unname(junc_genes[names(junc)])
  return(junc)
}
