# Rcppsva 0.1.0
Rcpp based Surrogate Variable Analysis and Rewrite by RcppEigen and RcppGSL.

# Rcppsva 0.1.2
* Update added original version of bump hunting for DMR(differential methylation region) by non-parametric methods(Permutations && Bootstrap) 
* Update in `helper.R` # 392 | 
> sum(abs(sum(ix[,3])) >= abs(colSums(ix[,-(1:7)])))/ L -> sum(abs(sum(ix[,3])) <= abs(colSums(ix[,-(1:7)])))/ L