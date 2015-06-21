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

