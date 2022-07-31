#' @export
#' @importFrom BiocGenerics unlist sort order
#' @importFrom S4Vectors split
get_first_exon <- function(gtf, label_only = TRUE){
  gtf_copy <- gtf
  gtf <- gtf[gtf$type=="exon"]
  gtf <- sort(gtf)
  minus <- gtf[strand(gtf)=="-"]
  minus <- unlist(last_element(split(minus, minus$transcript_id)))
  plus <- gtf[strand(gtf)=="+"]
  plus <- unlist(first_element(split(plus, plus$transcript_id)))
  first_exon <- c(plus, minus)
  if(label_only){
    gtf$first_exon <- paste0(gtf, "_", gtf$transcript_id) %in% paste0(first_exon, "_", first_exon$transcript_id)
    gtf_copy <- gtf_copy[gtf_copy$type != "exon"]
    if(length(gtf_copy) > 0){
      gtf_copy$first_exon <- FALSE
    }
    gtf <- c(gtf_copy, gtf)
    return(gtf)
  }
  first_exon <- first_exon[order(first_exon$transcript_id)]
  return(first_exon)
}

#' @export
#' @importFrom BiocGenerics unlist sort order
#' @importFrom S4Vectors split
get_last_exon <- function(gtf, label_only = TRUE){
  gtf_copy <- gtf
  gtf <- gtf[gtf$type=="exon"]
  minus <- gtf[strand(gtf)=="-"]
  minus <- unlist(first_element(split(minus, minus$transcript_id)))
  plus <- gtf[strand(gtf)=="+"]
  plus <- unlist(last_element(split(plus, plus$transcript_id)))
  last_exon <- c(plus, minus)
  if(label_only){
    gtf$last_exon <- paste0(gtf, "_", gtf$transcript_id) %in% paste0(last_exon, "_", last_exon$transcript_id)
    gtf_copy <- gtf_copy[gtf_copy$type != "exon"]
    gtf_copy$last_exon <- FALSE
    gtf <- c(gtf_copy, gtf)
    return(gtf)
  }
  last_exon <- last_exon[order(last_exon$transcript_id)]
  return(last_exon)
}

#' @export
get_single_exon <- function(x){
  if(class(x)=="character"){
    ext <- gsub(".*\\.","",basename(x))
    if(ext=="gtf"){
      x <- rtracklayer::import(x)
    }else{
      stop("only accept GTF format")
    }
  }
  exon <- x[x$type == "exon"]
  is_single <- table(exon$transcript_id)==1
  single_exon <- names(is_single)[is_single]
  return(single_exon)
}

#' @export
get_exon_grp <- function(anno, exon_type = "first", ess_tol, ees_tol){
  if(exon_type == "first"){
    exon <- get_first_exon(anno, label_only = FALSE)
  }
  if(exon_type == "last"){
    exon <- get_last_exon(anno, label_only = FALSE)
  }
  exon <- merge_exon(exon, ess_tol = ess_tol, ees_tol = ees_tol)
  exon <- exon[base::order(exon$transcript_id)]
  grp <- paste0(exon$gene_id, "_", as.character(exon))
  # ids <- paste0(exon$transcript_id, "_", exon$gene_id)
  members <- CharacterList(split(exon$transcript_id, grp))
  fe_grp <- rep(names(members), lengths(members))
  names(fe_grp) <- unlist(members, use.names = FALSE)
  fe_grp
}


#' @export
get_exon_count <- function(counts, gtf, sample_pair, ess_tol, ees_tol, exon_type){
  if(exon_type == "first"){
    grp <- get_exon_grp(gtf, exon_type = "first", ess_tol = 75, ees_tol = 0)
  }
  if(exon_type == "last"){
    grp <- get_exon_grp(gtf, exon_type = "last", ess_tol = 0, ees_tol = 75)
  }
  grp <- grp[rownames(counts)]
  members <- CharacterList(split(names(grp), grp))
  counts <- aggregate(counts[,unlist(sample_pair)], FUN = sum, by = list(grp))
  colnames(counts)[1] <- "coordinate"
  counts$members <- members[counts$coordinate]
  counts$symbols <- unlist(first_element(strsplit(counts[[1]], "_")), use.names = FALSE)
  counts
}
