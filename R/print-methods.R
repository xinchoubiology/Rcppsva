#' out.R is used to print/plot data object in Rcppsva to validation
#' output 450k microarray data frame with ordered positions to BED file
#' 
#' @rdname print
#' @param x m x n matrix. Statitical test p-value and each row's name is its probes name
#' @param ... bed BED file name
#'            db  Database of probes
#' @return BED data.frame
#' @method print BED
#' @examples 
#' require(IlluminaHumanMethylation450kanno.ilmn12.hg19)
#' data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
#' @export
print.BED <- function(x, bed = NULL, db = IlluminaHumanMethylation450kanno.ilmn12.hg19, ...){
  dat.m <- x$table
  Location <- data.frame(db@data$Locations)
  Probes <- rownames(Location)[rownames(Location) %in% rownames(dat.m)]
  BED <- cbind(Location[Probes, 1:2])
  BED <- cbind(BED, dat.m[Probes, 4:3])
  colnames(BED) <- c("CHR", "start", "pvalue", "F")
  
  BED <- GenomicRanges::GRanges(seqnames = BED$CHR, ranges = IRanges::IRanges(start = BED$start, width = 2), p = BED$pvalue, F = BED$F)
  seqlevels(BED) <- sort(seqlevels(BED))
  BED <- as.data.frame(sort(BED))[,c(1:3,6,7)]
  colnames(BED) <- c("chrom", "start", "end", "pvalue", "F")
  if(!is.null(bed)){
    write.table(as.matrix(BED), file = bed, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  print(head(BED, 5))
  cat(sprintf("... %d rows are omitted", nrow(BED) - 5))
}
