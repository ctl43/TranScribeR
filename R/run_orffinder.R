#' @export
run_orffinder <- function(input, output, tmp_dir = NULL){
  # Renaming sequence name as orffinder do not allow special character
  raw_fa <- read_fasta(input)
  renamed <- raw_fa
  names(renamed) <- seq_along(raw_fa)

  # Creating temp folder for the temp output
  if(is.null(tmp_dir)){
    tempdir()
    tmp_dir <- tempfile()
    dir.create(tmp_dir)
  }
  tmp_out_1 <- file.path(tmp_dir, paste0("renamed_", basename(input)))
  write_fasta(renamed, tmp_out_1)
  tmp_out_2 <- file.path(tmp_dir, paste0("temp_", basename(output)))

  # Run ORFfinder
  cmd <- paste(c("/home/ctlaw/tools/ORFfinder",
                 "-in", tmp_out_1,
                 "-s 0 -ml 9 -strand plus",
                 "-out", tmp_out_2), collapse = " ")
  system(cmd)

  # Convert back to the original names
  orff_out <-  read_fasta(tmp_out_2)
  names(orff_out) <- gsub(".*\\|", "", names(orff_out))
  tmp_1 <- strsplit(names(orff_out), "\\:")
  tmp_1 <- do.call(rbind, tmp_1)
  idx <- as.integer(gsub(".*_", "", tmp_1[, 1]))
  tmp_1[,1] <- gsub("_.*", "",tmp_1[, 1])
  names(orff_out) <- paste0(tmp_1[, 1], "|", names(raw_fa)[idx], "|", tmp_1[,2], "|", tmp_1[,3])
  names(orff_out) <- sub(" .*","",names(orff_out))
  write_fasta(orff_out, output)
  system(paste("rm -r", tmp_dir, collapse = " "))
}
