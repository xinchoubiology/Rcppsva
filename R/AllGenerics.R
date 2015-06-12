##' @docType methods
##' @name get.bumps
##' @rdname get.bumps-methods
##' @title get.bumps method
##' @param object methylation expression matrix object
##' @param ... additional parameters
##' @return bumps object
##' @export
setGeneric("get.bumps", function(object, ...){
  standardGeneric("get.bumps")
})

setOldClass("print")

##' @docType methods
##' @name print
##' @rdname print-methods
##' @title print method
##' @param object methylation expression matrix object
##' @param ... additional parameters
##' @export
setGeneric("print", function(object, ...){
  standardGeneric("print")
})



