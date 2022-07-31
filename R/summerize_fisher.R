#' @export
mark_diff <- function(x, min_frac = 0.03, pq="q", cutoff = 0.01, frac_fold = 2){
  .mark_diff <- function(y, min_frac, pq, cutoff, frac_fold){
    # Marking the differential test result
    y$trend <- "Unchanged"
    y$trend[y$frac_fold < (1/frac_fold)] <- "Decreased"
    y$trend[y$frac_fold > frac_fold] <- "Increased"
    y$trend[y[[pq]] > cutoff] <- "Unchanged"
    y$trend[y$frac_1 <= min_frac & y$frac_2 <= min_frac] <- "Unchanged"
    y$trend <- factor(y$trend, levels = c("Increased", "Decreased", "Unchanged"))
    y
  }
  lapply(x, .mark_diff, min_frac = min_frac, pq=pq, cutoff = cutoff, frac_fold = frac_fold)
}

#' @export
count_trend <- function(isoform_test, select = c("symbols", "transcript_id")){
  trends <- setDT(data.frame(sapply(isoform_test, "[[", "trend")))
  # For generating summary from a differential test
  trend_count <- data.frame(unchnaged = rowSums(trends == "Unchanged"),
                            increased = rowSums(trends == "Increased"),
                            decreased = rowSums(trends == "Decreased"))
  unique_in <- trend_count$increased
  unique_in[trend_count$decreased > 0] <- 0
  unique_de <- trend_count$decreased
  unique_de[trend_count$increased > 0] <- 0
  n_trend <- cbind(trend_count, uni_increased = unique_in, uni_decreased = unique_de)
  iso_trends <- cbind(trends, n_trend)
  iso_trends <- cbind(isoform_test[[1]][,select], iso_trends)
  iso_trends <- DataFrame(iso_trends)
  return(iso_trends)
}
