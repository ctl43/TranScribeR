#' @export
.fisher <- function(mat, selected, alt = NULL, grp = NULL, extra_info = NULL){
  # If grp is not NULL, this function will combine number in each column by group
  # Otherwise, it will compare the selected column and the alternative column by fisher.exact test.
  original_mat <- mat
  mat <- data.frame(mat)
  if(!is.null(grp)){
    grp <- mat[, grp]
    mat <- mat[, selected]
    original <- grp
    grp <- factor(grp)
    order <- order(grp)
    grp <- grp[order]
    mat <- mat[order,]
    mat <- setDT(mat)
    mat_copy <- mat
    colnames(mat) <- c("a","b")
    mat$grp <- grp
    grp_sum <- mat[, list(a = sum(a), b = sum(b)), by = grp]
    n <- table(grp)
    grp_sum <- do.call(cbind.data.frame, lapply(grp_sum[,-1], rep, times = n))
    alt <- data.table(grp_sum) - mat_copy
    combined <- cbind(mat_copy, alt)
  }else{
    if(!is.null(alt)){
      selected <- c(selected, alt)
      combined <- mat[,selected]
    }
  }

  combined <- t(combined)
  stat_test <- apply(combined, 2, function(x)fisher.test(matrix(x,2, byrow = TRUE)))
  p <- sapply(stat_test, "[[", "p.value")
  combined <- t(combined)
  colnames(combined) <- c("isoform_1", "isoform_2", "alt_1", "alt_2")
  out <- cbind.data.frame(combined, p)
  out <- DataFrame(out)
  out$frac_1 <- out$isoform_1/(out$isoform_1 + out$alt_1)
  out$frac_2 <- out$isoform_2/(out$isoform_2 + out$alt_2)
  if(!is.null(grp)){
    out <- out[order(order),]
  }
  out$frac_fold <- out$frac_2/out$frac_1
  if(!is.null(extra_info)){
    out <- cbind(out, original_mat[,extra_info])
  }
  return(out)
}

#' @export
#' @importFrom BiocParallel SerialParam bpmapply bplapply
fisher <- function(mat, sample_pair, selected, alt = NULL, grp = NULL, extra_info = NULL, BPPARAM = SerialParam()){
  if(class(sample_pair) == "character" & length(sample_pair) == 2){
    sample_pair <- list(sample_pair)
    retrun(bplapply(sample_pair, .fisher, mat = mat, grp = grp, BPPARAM = SerialParam())[[1]])
  }
  if(class(sample_pair) == "list"){
    if(class(alt)=="list"){
      return(bpmapply(function(x,y).fisher(x,y, mat = mat, grp = grp, extra_info = extra_info),
                      x = sample_pair, y = alt, BPPARAM = BPPARAM, SIMPLIFY = FALSE))
    }
    return(bplapply(sample_pair, .fisher, mat = mat, grp = grp, BPPARAM = BPPARAM, extra_info = extra_info))
  }
}

