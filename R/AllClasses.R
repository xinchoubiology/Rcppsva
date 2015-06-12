setOldClass("print")

##' Class "bumps"
##' The class contains information of methyaltion differential related region
##' 
##' @name bumps-class
##' @docType class
##' @slot table statistic table of CpG regions
##' @slot null list of null distribution (H0: no differential methylation)
##' @slot algorithm "get.bumps" algorithm's argumrnts
##' @exportClass bumps
##' @keywords classes 
setClass("bumps",
         representation = representation(bumps     = "list",
                                         algorithm = "list")
)


                        
                        