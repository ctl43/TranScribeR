# predict_poison <- function(x, bed, paried_by = "cds"){
#   # srsf <- orf[grep("SRSF", mcols(orf)$transcript_ids)]
#   # srsf <- srsf[mcols(srsf)$main_orf]
#   # mcols(srsf)$gene <- unname(paste(CharacterList(first_element(strsplit(mcols(srsf)$transcript_ids, "_"), invert = TRUE)), collapse = "_"))
#   # nmd <- srsf[mcols(srsf)$is_nmd]
#   # nmd <- split(mcols(nmd)$transcript_ids, mcols(nmd)$gene)
#   # nmd_srsf3 <- nmd[["SRSF3"]]
#   #
#   #
#   #
#   # bed <- import_bed("flair_out/collapsed/flank_100_only/assigned_collapsed.isoforms.bed")
#
#   gtf <- rtracklayer::import("/home/ctlaw/reference/Homo_sapiens/UCSC_hg19/Annotation/gencode.v38lift37.annotation.gtf")
#   table(gtf$type)
#   gtf <- rtracklayer::import("flair_out/collapsed/flank_100_only/")
#   gtf$uid <- paste0(gtf$transcript_id, "_", gtf$gene_id)
#   gtf <- gtf[gtf$type=="exon"]
#   srsf3_gtf <- gtf[gtf$gene_id=="SRSF3"]
#   srsf3_exon <- split(as.character(srsf3_gtf), srsf3_gtf$uid)
#   srsf3_exon <- CharacterList(srsf3_exon)
#
#   # srsf3_bed <- bed[mcols(srsf[mcols(srsf)$gene=="SRSF3"])$transcript_ids]
#   is_nmd <- names(srsf3_exon)%in%nmd_srsf3
#
#   srsf3_nmd <- srsf3_exon[is_nmd]
#   srsf3_non_nmd <- srsf3_exon[!is_nmd]
#
#
#   srsf3_bed$non_nmd <- srsf3_bed$non_nmd[srsf3_bed$non_nmd%in%srsf3_bed$nmd]
#   x <- srsf3_nmd["127:1511|4e261e3b-539c-450b-8bd2-0eb6c6ac519c-1_SRSF3"]
#   non_nmd <- srsf3_non_nmd#[["119:2049|69056042-9354-48a3-950a-e15d4af9be23_SRSF3"]]
#   # Rmb to use coding region for comparison only!!!!!!!!
#
#   get_poison(srsf3_nmd, srsf3_non_nmd)
#   get_poison <- function(nmd, non_nmd){
#     compare_exon <- function(x, non_nmd){
#       y <- CharacterList(x)
#       is_paired <- sum(!y%in%non_nmd)==1
#       if(any(is_paired)){
#         return(non_nmd[is_paired][[1]])
#       }else{
#         return(x)
#       }
#     }
#     paired <- CharacterList(lapply(srsf3_nmd, compare_poison, non_nmd = non_nmd))
#     return(paired)
#   }
#
#
#
#   poison <- srsf3_nmd[!srsf3_nmd%in%paired]
#
# }
