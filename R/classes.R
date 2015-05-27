setClass("bumps",
         representation(table     = "data.frame",
                        coef      = "numeric",
                        fitted    = "numeric",
                        pvalue    = "numeric",
                        null      = "list",
                        algorithm = "list"))

setMethod("print", signature(x = "bumps"),
          function(x, ...){
            cat(sprintf("'bumps' object with %d bumps\n", ncol(x@table)))
            cat(sprintf("algorithms : ...... \n"))
            for(params in names(x@algorithm)){
              cat(sprintf("\t%s = %s \n", params, x@algorithm[[params]]))
            }
          })
                        
                        