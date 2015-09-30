#' IMPLEMENTATOIN OF SMMB TEST OPTION OF GWJ
#'
#' Computes tests statistics for weights calculated according the smmb option
#' Internal function called by GWJsieve function \emph{not} gonna
#' export it.
#'
#' @param weight.obj A matrix.
#' @param vaccine.indices indices into \code{weight} for vaccine recipients.
#' @param placebo.indices indices into \code{weight} for placebo recipients.
#' @param ... additional arguments
smmb.test <- function(weight.obj,vaccine.indices , placebo.indices,...){
  
  #NOTE: CASTING INTO THE MATRIX FORMAT IS DONE SO THAT IF WE ONLY HAVE ONE ROW IT REMAINS 1XN MATRIX AND NOT VECTOR OF SIZE N
  v.count <- apply(matrix(weight.obj$weights.count[vaccine.indices,],nrow=length(vaccine.indices)),MARGIN=2,FUN=sum)
  p.count <- apply(matrix(weight.obj$weights.count[placebo.indices,],nrow=length(placebo.indices)),MARGIN=2,FUN=sum)
 
  ###THIS PART PERHAPS CAN BE REMOVED. IT JUST ENSURES NA'S SHOW IN PLACES WHERE THE COLUMN WAS 0 TO MATCH THE OLD ALGORITHM, WE COULD HAVE ZERO THERE
  isValid <- rep(1,ncol(weight.obj$weights))
  isValid[which((apply(weight.obj$weights[vaccine.indices,],MARGIN=2,FUN=sd) == 0) & (apply(weight.obj$weights[placebo.indices,],MARGIN=2,FUN=sd) == 0) ) ] <- NA
  ###
  
  phat.v <- apply(matrix(weight.obj$weights[vaccine.indices,],nrow=length(vaccine.indices)), MARGIN=2,FUN=sum)/v.count
  phat.p <- apply(matrix(weight.obj$weights[placebo.indices,],nrow=length(placebo.indices)), MARGIN=2,FUN=sum)/p.count
  return (list(test.stats = (phat.v - phat.p)*isValid ))
}