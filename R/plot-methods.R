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
  mat.m$value <- cut(mat.m$value,breaks=c(-1.01,-abs(cutoff), 0, abs(cutoff), 1.01),include.lowest=FALSE,label=c("TRUE - ","FALSE - ","FALSE + ","TRUE + "))
  mat.m$row   <- factor(mat.m$row, levels=rev(unique(as.character(mat.m$variable))))
  ggcor       <- ggplot(mat.m, aes(row, variable)) +
                    geom_tile(aes(fill=value),colour="black")  +
                    scale_fill_brewer(palette = "RdYlGn",name="Correlation") + 
                    theme(axis.text.x=element_text(angle=-90), 
                          axis.title=element_blank()) +
                    theme(panel.background=element_blank(),
                          panel.grid.minor=element_blank(),
                          panel.grid.major=element_blank())
  ggcor
}

#' plot hierachical tree for clustering data by ggplot
#' @rdname plot
#' @param x dendro class object
#' @param ... ignored argument
#' @param rotate if TRUE, rotate tree by 90 degree; Default FALSE
#' @param labels if TRUE, show segment labels
#' @method plot dendro
#' @export
plot.dendro <- function(x, ..., rotate = FALSE, labels = TRUE){
  angle <- ifelse(rotate, 0, 90)
  hjust <- ifelse(rotate, 0, 1)
  p <- ggplot() + 
       geom_segment(data = x$segments, aes_string(x = "x", y = "y", xend = "xend", yend = "yend"))
  if(labels){
    p <- p + scale_x_discrete(labels = x$labels$label)
  }
  if(rotate){
    p <- p + coord_flip()
    p <- p + scale_y_continuous()
  }else{
    p <- p + scale_y_continuous()
  }
  p <- p + theme_dendro()
  p <- p + theme(axis.text.x = element_text(angle = angle, hjust = 1))
  p <- p + theme(axis.text.y = element_text(angle = angle, hjust = 1))
  p
}


#' create complete blank theme in ggplot for dendro class object
#' @title theme_dendro
#' @export
theme_dendro <- function(){
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x     = element_text(colour = NA),
        axis.title.y     = element_blank(),
        axis.text.x      = element_blank(),
        axis.text.y      = element_blank(),
        axis.line        = element_blank(),
        axis.ticks       = element_blank())
}
