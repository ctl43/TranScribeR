# TranScribeR
This package provide tools for analysing unannotated transcriptome.

## Comparing assembled annotation and reference
This method only compares the junction between isoform. Regardless the end of a transcript, if the splice junctions of a transcript completely match with the splice junctions of a reference transcript, it is classified as full-splice matches (FSM); if they only partially matches the junctions of a reference transcript, it is classified as incomplete-splice matches (ISM); if they have a novel connection between annotated junctions, it is classified as novel in catalog; if they have an unannotated splice junction, it is classified as novel not in catalog. This classification is reference to SQANTI (Tardaguila, et al., 2018).

### Example
```r
junctions <- get_junctions(gtf) # newly assmebled transcriptome
ref_junctions <- get_junctions(gencode, gene_id = "gene_name") # reference transcriptome
out <- compare_junction(junctions, ref_junctions, identifier = "genes", BPPARAM = MulticoreParam(workers = 2L))
```

## Predicitng ORF
### Getting ORF (`predict_orf` and `get_main_orf`)
`predict_orf` predicts all open reading frames (ORF) from a sequences based on filters, for examples, 
start codon, stop 
codon and the minimum length of proteins, etc.
If reference protein sequence is supplied (e.g. uniprot), the protein from the predicted ORF will be renamed to the reference protein when there is any matched amino acid sequence in the reference annotation.
`get_main_orf` integrates the result from hmmscan/hmmsearch and blastp to determine the main ORF. 
Additionally, it predicts whether the mRNA will be degraded by nonsense-mediated mRNA decay (NMD) if 
it fulfill certain requirements (See Section, Predicting NMD).

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

## Predicting NMD
###Example
```r
nmd_predicted <- predict_nmd(tx_orf = morf$tx_coord, anno = anno)
morf$nmd <- nmd_predicted
```

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

## Format manupulation
### Example
```r
bed12 <- import_bed12(bed12_file)
gtf <- bed122gtf(bed12)
bed12 <- gtf2bed12(gtf)
export_bed12(bed12, out = path)
```

## Coordinate conversion
### Example
```r
# Converting from genomic coordinate to transcript coordinate
# For example, converting the coordinate of a start and stop codon in genome to the coordinate in transcript
tx_coord <- genome2tx(anno = protein_coding_annotation)

# Converting from transcript coordinate to genomic coordinate
genomic_coord <- tx2genome(tx_coord, anno = bed12)
```
