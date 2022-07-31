#' @export
#' @importFrom reshape2 melt
#' @importFrom igraph graph_from_data_frame components
merge_exon <- function(exon, ess_tol = 150, ees_tol = 10){
  ess <- resize(exon, width = 1)
  ees <- resize(exon, width = 1, fix = "end")

  ess_ol <- findOverlaps(ess, ess + ess_tol)
  ess_ig <- graph_from_data_frame(as.data.frame(ess_ol))
  ess_m <- components(ess_ig)$membership
  ess_m <- split(as.integer(names(ess_m)), ess_m)
  ess_m <- melt(ess_m)
  ess_m <- ess_m[order(ess_m[,1]),]

  ees_ol <- findOverlaps(ees, ees + ees_tol)
  ees_ig <- graph_from_data_frame(as.data.frame(ees_ol))
  ees_m <- components(ees_ig)$membership
  ees_m <- split(as.integer(names(ees_m)), ees_m)
  ees_m <- melt(ees_m)
  ees_m <- ees_m[order(ees_m[,1]),]

  grp <- cbind(ess_m, ees_m[,2])
  colnames(grp) <- c("idx", "ess_grp", "ees_grp")
  combined_grp <- paste0(grp$ess_grp, "_", grp$ees_grp)
  new_ess <- resize(range(split(ess[grp[[1]]], combined_grp)), width = 1)
  new_ees <- resize(range(split(ees[grp[[1]]], combined_grp)), width = 1, fix = "end")
  new_ess <- unlist(new_ess)
  new_ees <- unlist(new_ees)
  member <- IntegerList(split(grp[[1]], combined_grp))

  member <- member[unique(combined_grp)]
  new_ess <- new_ess[unique(combined_grp)]
  new_ees <- new_ees[unique(combined_grp)]

  new_ess <- rep(new_ess, lengths(member))
  new_ees <- rep(new_ees, lengths(member))
  new_ess <- new_ess[order(unlist(member))]
  new_ees <- new_ees[order(unlist(member))]

  start(exon[strand(exon)=="+"]) <- start(new_ess[strand(exon)=="+"])
  end(exon[strand(exon)=="-"]) <- start(new_ess[strand(exon)=="-"])
  end(exon[strand(exon)=="+"]) <- start(new_ees[strand(exon)=="+"])
  start(exon[strand(exon)=="-"]) <- start(new_ees[strand(exon)=="-"])

  return(unname(exon))
}
