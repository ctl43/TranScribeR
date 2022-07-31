# get_covered_exon <- function(target, bed, gtf){
#   bed <- "/home/ctlaw/dicky/analysis/nanopore/combined_hypoxia/flair_out/collapsed/assigned_collapsed.isoforms.bed"
#   library(TranScribeR)
#   if(class(bed)=="character"){
#     bed <- import_bed(bed)
#     anno <- bed122gtf(bed)
#   }else{
#     if(class(bed)!="GRanges"){
#       stop()
#     }else{
#       anno <- bed122gtf(bed)
#     }
#   }
#
#   if(class(gtf)=="character"){
#     anno <- import(gtf)
#   }else{
#     if(class(gtf)!="GRanges"){
#       stop()
#     }else{
#       anno <- gtf
#     }
#   }
#
#   exon <- anno[anno$type=="exon"]
#   exon <- split(exon, paste0(exon$transcript_id, "_", exon$gene_id))
#   exon <- exon[target$transcript_id]
#   cds_exon <- mcmapply(restrict, exon, start = start(target), end = end(target), mc.cores = 10L) #should be optimised
#   cds_exon <- GRangesList(cds_exon)
#   cds_exon
# }
