#' @export
check_novelty <- function(q, r){
  collected <- NULL
  status <- ""
  q_len <- length(q)
  r_len <- lengths(r)
  for(i in seq_along(r)){
    current <- r[[i]]
    which_ol <- which(current == q[1])
    current_len <- r_len[i]
    if(length(which_ol) == 0){
      next
    }else{
      if((which_ol + q_len - 1) > current_len){
        next
      }
      selected <- current[which_ol:(which_ol + q_len - 1)]
      all_match <- all(selected == q)
      if(all_match){
        if(current_len == q_len){
          if(status != "FSM"){
            collected <- c()
            status <- "FSM"
          }
          if(status == "FSM"){
            collected <- c(collected, c(names(r)[i]))
          }
        }else{
          if(status == "FSM"){
            next
          }
          status <- "ISM"
          if(status == "ISM"){
            collected <- c(collected, c(names(r)[i]))
          }
        }
      }else{
        next
      }
    }
  }
  if(status == ""){
    total_ss <- unlist(r)
    if(all(as.logical(BiocGenerics::match(q, total_ss, nomatch = 0)))){
      return(c("", "NIC"))
    }else{
      return(c("", "NNC"))
    }
  }else{
    return(c(paste(collected, collapse = ","), status))
  }
}

#' @export
compare_junction <- function(q, ref_junctions, BPPARAM = SerialParam(), identifier = "gene_name"){
  # This limited the comparison within isoforms that are assigned to the same gene.
  q_gene <- mcols(q)[[identifier]]
  list_grp <- get_list_grp(q)
  r_gene <- mcols(ref_junctions)[[identifier]]
  q <- sort(q)
  ref_junctions <- sort(ref_junctions)
  out <- bpmapply(function(a,b){
    r <- ref_junctions[r_gene == b]
    check_novelty(a, r)
  }, a = q, b = q_gene, BPPARAM = BPPARAM)
  out <- t(out)
  out <- cbind(names(q), out)
  out <- DataFrame(out)
  colnames(out) <- c("query", "related_ref", "type")
  out[[2]] <- CharacterList(strsplit(out[[2]], ","))
  rownames(out) <- NULL
  return(out)
}
