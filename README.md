# TranScribeR
This package provide useful tools and wrapper for transcriptomic analysis.

## Predicitng ORF
### Getting ORF (`predict_orf` and `get_main_orf`)
`predict_orf` predicts all open reading frames (ORF) from a sequences based on filters, for examples, 
start codon, stop 
codon and the minimum length of proteins, etc.
If reference protein sequence is supplied (e.g. uniprot), the protein from the predicted ORF will be renamed to the reference protein when there is any matched amino acid sequence in the reference annotation.
`get_main_orf` integrates the result from hmmscan/hmmsearch and blastp to determine the main ORF. 
Additionally, it predicts whether the mRNA will be degraded by nonsense-mediated mRNA decay (NMD) if 
it fulfill certain requirements (See Section, Predicting transcription outcome).

### Example
```r
# Predicting ORF and assigning protein ID
predicted <- predict_orf(fa, seq = NULL, start_codon = "ATG", starts = NULL, 
		stop_codon = c("TAG", "TAA", "TGA"), min_aa = 30, ref_prot_fa = NULL, 
		assign_prot_id = TRUE, with_stop = FALSE) 

# Predicting the main ORF out of many predicted ORF
orf <- get_main_orf(predicted_orf = predicted,
                    anno = isoform_annotation_bed12,
                    pfam_out = NULL, 
                    blastp_out = NULL,
                    ref_start_codon = NULL, 
                    min_aa = 50, 
                    label_main_only = TRUE, 
                    predict_nmd = TRUE)
morf <- orf[orf$main_orf,]

# Appending those non-coding mRNA
morf <- append_noncoding(morf, fa = fa)
```

## Predicting regulatory transcription outcome

## Predicting upstream ORF (uORF) (`get_uorf`)
It predicts uORF after getting the main ORF from above command or other methods. In addition, users can supply reference start site retrieved by ribo-seq or online annotations to constraint the uORF prediction.

### Example
```r
uorf <- get_uorf(morf$tx_coord, morf$phase, fa, anno = bed12, 
         ref_start, 
         start_codons = c("ATG", "TTG", "ACG", "CTG", "GTG"), 
         stop_codons = c("TAG", "TAA","TGA"), 
         min_aa = 2)
```

