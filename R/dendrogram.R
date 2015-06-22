## create dendrogram df object for hierachical clustering structure building

#' extract cluster data frame from a model into a list of data frame
#' @title dendro.data
#' @param model object of \code{\link{stats}{hclust}}
#' @param ... ignored
#' @aliases dendro.data.default
#' @rdname dendro.data
#' @export
dendro.data <- function(model, ...){
  UseMethod("dendro.data", model)
}

#' default dendro data convert 
#' dendro.data method for \code{default} object
#' @rdname dendro.data
#' @param model model object 
#' @return NULL
#' @method dendro.data default
#' @export
dendro.data.default <- function(model, ...){
  x <- class(model)
  stop(paste("hclust class needed by dendro.data method but not", x))
  return(NULL)
}

#' dendro data convert for \code{\link{stats}{hclust}} object
#' @rdname dendro.data
#' @param model model object; hc
#' @return ggplot of clustering hierachical tree
#' @method dendro.data hclust
#' @export
dendro.data.hclust <- function(model, type = c("rectangle", "triangle"), ...){
  type <- match.arg(type)
  dendro <- as.dendrogram(model)
  res <- dendrogram(dendro, type = type)
  return(res)
}


#' extracting data frame from hclust object for plotting dendrgram
#' @title dendrogram
#' @param x object of hclust. Derived from \code{\link{stats}{dendrogram}}
#' @param type type of plot, c("rectangle", "triangle"); Default rectangle
#' @param ... ignore
#' @return list
#'         segments data.frame
#'         labels data.frame
#'         class
#' @export
dendrogram <- function(x, type = c("rectangle", "triangle"), ...){
  x1 <- 1
  x2 <- .memberDend(x)
  type <- match.arg(type)
  dendronode <- function(x1, x2, subtree, ddsegments = NULL, ddlabels = NULL){
    inner <- !is.leaf(subtree) && x1 != x2
    bx <- .memberLimit(x1, x2, subtree)
    x.top <- bx$x
    y.top <- attr(subtree, "height")
    if(is.leaf(subtree)){
      ddlabels <- cbind(ddlabels, data.frame(x = x.top, y = 0, label = attr(subtree, "label")))
    }else if(inner){
      for(k in 1:length(subtree)){
        child <- subtree[[k]]
        y.bot <- attr(child, "height")
        x.bot <- bx$limit[k:(k+1)] + .midDend(child)
        if(type == "rectangle"){
          ddsegments <- cbind(ddsegments, data.frame(x = x.top, y = y.top, xend = x.bot, yend = y.bot))
        }else{
          ddsegments <- cbind(ddsegments, data.frame(x = x.top, y = y.top, xend = x.bot, yend = y.top))
          ddsegments <- cbind(ddsegments, data.frame(x = x.bot, y = y.top, xend = x.bot, yend = y.bot))
        }
        res <- dendronode(x1 = bx$limit[k], x2 = bx$limit[k+1], child, ddsegments = ddsegments, ddlabels = ddlabels)
        ddsegments <- res$segments
        ddlabels <- res$labels
      }
    }
    return(list(segments = ddsegments, labels = ddlabels))
  }
  res <- dendronode(x1, x2, x, ddsegments = NULL, ddlabels = NULL)
  names(res$segments) <- c("x", "y", "xend", "yend")
  names(res$labels)   <- c("x", "y", "label")
  res
}


#' assign position of fork and two limit of tree to members
#' @title .memberLimit 
#' @param x1 leftmost node 
#' @param x2 rightmost node
#' @param subtree subtree of dendrogram
#' @return list
#'         x
#'         limit c(tree1, tree2, tree3)
.memberLimit <- function(x1, x2, subtree){
  if(!is.leaf(subtree) && x1 != x2){
   k <- length(subtree)
   limit <- integer(k)
   xx <- x1
   for(i in 1:k){
     m <- .memberDend(substree[[i]])
     xx <- xx + m
     limit[i] <- xx
   }
   limit <- c(x1, limit)
  }else{
   limit <- c(x1, x2) 
  }
  mid <- attr(subtree, "midpoint")
  if(!is.numeric(mid)){
    x <- x1 + (if(!is.leaf(subtree) && x1 != x2)
                 mid
               else
                  0)
  }else{
    x <- mean(c(x1, x2))
  }
  
  return(list(x = x, limit = limit))
}

#' return member number in subtree x
#' @param subtree
#' @return number in tree object x
.memberDend <- function(subtree){
  r <- attr(subtree, "members")
  if(is.null(r))
    r <- 1L
  r
}

#' return midpoint of subtree
#' @param subtree
#' @return midpoint value
.midDend <- function(subtree){
  mp <- attr(subtree, "midpoint")
  if(is.null(mp)){
    mp <- 0
  }
  mp
}
