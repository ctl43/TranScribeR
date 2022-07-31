# TranScribeR
This package provide useful tools and wrapper for transcriptomic analysis.

## Predicitng ORF
### Running ORFfinder (`run_orffinder`)
A wrapper function for calling ORFfinder, which is provided by NCBI, from the sysetm.

### Getting ORF (`get_orf`)
Result from ORFfinder was imported and get the longest predicted ORF as the final ORF if there are more than one ORF in a trascript.
When results from hmmsearch/hmmscan and blastp are available, it will also integrate these information in predicting the final ORF.
If reference protein sequence is supplied (e.g. uniprot), the protein from the predicted ORF will be renamed to the reference protein when there is any matched amino acid sequence in the reference annotation.
The program does not do any complicated algorithm, (e.g. TransDecoder) to integrate the result from hmmscan/hmmsearch and blastp.
It just simply marks the ORF as potetinal candidates when protein encoded from the predicted ORF has any domain/similar protein suggested by either hmmscan/hmmsearch or blastp.

### Example
```r
run_orffinder(input = "seq.fa", output = "orf.fa")
orf <- get_orf(orffinder_out = "orf.fa", fa = "seq.fa", pfam_out = "hmmsearch.domtblout", 
               blastp_out = "blastp.outfmt6", ref_prot_fa = "uniprot.fasta", min_len = 50)
write_fasta(as.character(orf[mcols(orf)$main_orf]), "main_orf.fa")
```
