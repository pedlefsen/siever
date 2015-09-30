#' ALTERNATE IMPLEMENTATION OF GWJ Kullbak-Leibler stat driver
#'
#' Blah.  Blah.  Don't know what to say (Ted).  Looks like a very abstract implementation
#' of the algorithm.  Seems to depend purely on a weight matrix.  Someone else must write this.
#' Not sure if this is an internal function.  It looks like one.  Temporarily \emph{not} gonna
#' export it.
#'
#' @param weight.obj A matrix.
#' @param vaccine.indices indices into \code{weight} for vaccine recipients.
#' @param placebo.indices indices into \code{weight} for placebo recipients.
#' @param ... additional arguments
klStat.alt <- function(weight.obj,vaccine.indices,placebo.indices, ...){
  acceptable.chars <- weight.obj$opt$acceptable.chars
  symetric.stat = weight.obj$opt$stat.type == "SKL"
  data.mtx <- weight.obj$weights
  
  vaccine.props <- apply(data.mtx[vaccine.indices,],MARGIN=2,FUN=function(col){tab <- table(factor(col,levels=acceptable.chars)); return(tab/sum(tab))})
  placebo.props <- apply(data.mtx[-vaccine.indices,],MARGIN=2,FUN=function(col){tab <- table(factor(col,levels=acceptable.chars)); return(tab/sum(tab))})
  
  num.vaccine.seqs <- length(vaccine.indices)
  num.placebo.seqs <- nrow(data.mtx) - num.vaccine.seqs
  
  seq.length <- ncol(data.mtx)
  
  discrepancy <- function(a, b){ a * log(a/b) }
  if(symetric.stat){
    discrepancy <- function(a, b){ a * log(a/b) + b * log(b/a) }
  }
   
  kl.scores <- matrix(0,nrow= length(acceptable.chars), ncol = seq.length, dimnames=list(acceptable.chars,NULL))
  props.product.notZero <- vaccine.props * placebo.props != 0
  kl.scores[props.product.notZero] <- discrepancy(vaccine.props, placebo.props)[props.product.notZero]
  kl.scores[!props.product.notZero] <- discrepancy(vaccine.props + 1/num.vaccine.seqs,placebo.props + 1/num.placebo.seqs)[!props.product.notZero]

  test.stats <- apply(kl.scores, MARGIN=2,FUN=sum)
  
  return(list(test.stats = test.stats))

}