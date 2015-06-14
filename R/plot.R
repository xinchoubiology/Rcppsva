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

#' plot correlation matrix by ggplot
#'
#' @rdname plot
#' @param x n x n correlation matrix
#' @param y NULL
#' @param cutoff correlation cutoff; Default = 0.8
#' @return ggplot object
#' @method plot corr
#' @importFrom reshape2 melt
#' @import ggplot2
#' @import scales
#' @import RColorBrewer
#' @examples
#' corrmat <- list(cor = matrix(c(-0.9,-0.6,0.3,1), nrow = 2))
#' colnames(corrmat$cor) <- rownames(corrmat$cor) <- c("X1","X2")
#' class(corrmat) <- "corr"
#' plot(corrmat)
#' @export
plot.corr <- function(x, y = NULL, cutoff = 0.8, ...){
  if(!is.null(y)){
    stop("y object in plot.corr must be NULL")
  }
  mat <- data.frame(row = rownames(x$cor), x$cor)
  rownames(mat) <- NULL
  mat.m <- melt(mat)
  mat.m$value<-cut(mat.m$value,breaks=c(-1.01,-abs(cutoff), 0, abs(cutoff), 1.01),include.lowest=FALSE,label=c("TRUE - ","FALSE - ","FALSE + ","TRUE + "))
  mat.m$row <- factor(mat.m$row, levels=rev(unique(as.character(mat.m$variable))))
  ggcor <- ggplot(mat.m, aes(row, variable)) +
            geom_tile(aes(fill=value),colour="black")  +
            scale_fill_brewer(palette = "RdYlGn",name="Correlation") + 
            theme(axis.text.x=element_text(angle=-90), 
                  axis.title=element_blank()) +
            theme(panel.background=element_blank(),
                  panel.grid.minor=element_blank(),
                  panel.grid.major=element_blank())
  ggcor
}