#' Acumulating function over observations
#' @param weight.obj A weight object
obs.accumulate <- function(weight.obj){
  seq.length <- ncol(weight.obj$weights)
  vIndex <- weight.obj$v.indices
  pIndex <- (1:nrow(weight.obj$weights))[-vIndex]
  weights.v <- NULL
  weights.p <- NULL
  
  weights.count.v <- NULL
  weights.count.p <- NULL
  
  # ACCUMULATES WEIGHTS ACROSS OBSERVERATIONS IN THE VACCINE OBSERVATION SEQUENCE SETS
  if(!is.null(weight.obj$opt$obs.vSeq.sets.list)) {
    #obs.num.vSeq <- length(unique(unlist(weight.obj$opt$obs.vSeq.sets.list)))
    num.vObs.set <- length(weight.obj$opt$obs.vSeq.sets.list)
    weights.v <-
      matrix(unlist(
        lapply(weight.obj$opt$obs.vSeq.sets.list,
               FUN=function(set) apply( matrix(weight.obj$weights[set,],ncol=seq.length),MARGIN=2, FUN= weight.obj$opt$within.obs.acc.fn))),
        nrow = num.vObs.set,
        byrow=T
      )
    
    
    vObs.set.names <- names(weight.obj$opt$obs.vSeq.sets.list)
    if(is.null(vObs.set.names)) vec.vSet.names <- paste("vObs.set", 1:num.vObs.set, sep="")
    rownames(weights.v) <- vec.vSet.names
    
    # SMMB WEIGHTS COUNT ACCUMULATION
    if(!is.null(weight.obj$weights.count)){
      weights.count.v <-
        matrix(unlist(
          lapply(weight.obj$opt$obs.vSeq.sets.list,
                 FUN=function(set) apply( matrix(weight.obj$weights.count[set,],ncol=seq.length),MARGIN=2, FUN= weight.obj$opt$within.obs.acc.fn))),
          nrow = num.vObs.set,
          byrow=T
        )
      rownames(weights.count.v) <- vec.vSet.names
    }
    
    vIndex <- NULL
    
  }
  
  # ACCUMULATES WEIGHTS ACROSS OBSERVERATIONS IN THE PLACEBO OBSERVATION SEQUENCE SETS
  if(!is.null(weight.obj$opt$obs.pSeq.sets.list)) {
    #obs.num.pSeq <- length(unique(unlist(weight.obj$opt$obs.pSeq.sets.list)))
    num.pObs.set <- length(weight.obj$opt$obs.pSeq.sets.list)
    weights.p <-
      matrix(
        unlist(lapply(weight.obj$opt$obs.pSeq.sets.list,
                      FUN=function(set) apply( matrix(weight.obj$weights[set,],ncol=seq.length), MARGIN=2, FUN= weight.obj$opt$within.obs.acc.fn))
        ),
        nrow = num.pObs.set,
        byrow=T
      )
    
    # SMMB WEIGHTS COUNT ACCUMULATION
    pObs.set.names <- names(weight.obj$opt$obs.pSeq.sets.list)
    if(is.null(pObs.set.names)) vec.pSet.names <- paste("pObs.set", 1:num.pObs.set, sep="")
    rownames(weights.p) <- vec.pSet.names
    
    if(!is.null(weight.obj$weights.count)){
      weights.count.p <-
        matrix(unlist(
          lapply(weight.obj$opt$obs.pSeq.sets.list,
                 FUN=function(set) apply( matrix(weight.obj$weights.count[set,],ncol=seq.length),MARGIN=2, FUN= weight.obj$opt$within.obs.acc.fn))),
          nrow = num.vObs.set,
          byrow=T
        )
      rownames(weights.count.p) <- vec.pSet.names
    }
    
    pIndex <- NULL
  }  
  weights <- rbind(weight.obj$weights[vIndex,], weights.v, weight.obj$weights[pIndex,], weights.p)
  v.indices <- 1:(length(vIndex) + max(0, nrow(weights.v)))
  
  weights.count <- rbind(weight.obj$weights.count[vIndex,], weights.count.v, weight.obj$weights.count[pIndex,],weights.count.p)
  return(list(weights = weights, weights.count = weights.count ,v.indices = v.indices))
  
}
