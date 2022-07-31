#' @export
#' @importFrom rtracklayer import export
mapping_trans2gene <- function(input, output, reference, gffcompare_path = "gffcompare", tmp_dir = NULL){
  if(is.null(output)){
    output <- file.path(getwd(), "assigned")
  }
  # tmp_out <- file.path(tmp_dir, paste0("temp_", basename(output)))
  cmd <- paste(c(gffcompare_path, "-R -r", reference, "-o", output, input), collapse = " ")
  system(cmd)
  anno_copy <- anno <- import(paste0(output, ".annotated.gtf"))
  anno <- anno[anno$type == "transcript"]
  anno[is.na(anno$gene_name)]$gene_name <- gsub("_", "",anno[is.na(anno$gene_name)]$xloc)
  tr2gene <- anno$gene_name
  names(tr2gene) <- anno$transcript_id
  anno_copy$gene_id  <- tr2gene[anno_copy$transcript_id]
  mcols(anno_copy) <- mcols(anno_copy)[c("source","type","score","phase","gene_id","transcript_id", "exon_number")]
  export(anno_copy, paste0(output, "_", basename(input)))
  paste("rm", "-r", tmp_dir, collapse = " ")
  system(paste0("rm ", output, ".stats"))
  system(paste0("rm ", output, ".annotated.gtf"))
  system(paste0("rm ", output, ".loci"))
  system(paste0("rm ", output, ".tracking"))
  system(paste0("rm ", output, ".", basename(input), ".refmap"))
  system(paste0("rm ", output, ".", basename(input), ".tmap"))
  return(tr2gene)
}
