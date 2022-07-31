blastp_makedb <- function(makeblastdb, ref_aa, title = "reference_aa", out){
  cmd <- paste(c(makeblastdb, "-in", ref_aa, 
          "-input_type fasta", 
          "-dbtype prot",
          "-title", title, 
          "-parse_seqids -out", out), 
        collapse = " ")
  system(cmd)
}

run_blastp <- function(blastp, aa_fa, out, threads = 1){
  cmd <- paste(c(blastp, "-query", aa_fa, "-db", blastp_db, 
                        "-max_target_seqs", 1,
                        "-outfmt 6 -evalue 1e-5", 
                        "-num_threads", threads, ">", out),
                      collapse = " ")
  system(cmd)
}

run_hmmsearch <- function(hmmsearch, aa_fa, out, pfam_hmm, threads = 1){
  cmd <- paste(c(hmmsearch, "--cpu", threads, "--domtblout", out, 
                 pfam_hmm, aa_fa),
               collapse = " ")
  system(cmd)
}
