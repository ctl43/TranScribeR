#' @export
get_canonical_isoform <- function(x, factor){
  count <- rowSums(as.matrix(x))
  ids <- rownames(x)
  dt <- data.table(ids, factor, count)
  dt <- dt[,.SD[which.max(count),], by = factor]
  DataFrame(canonical = dt[["ids"]],
            factor = dt[["factor"]],
            members = CharacterList(split(ids, factor)))
}
