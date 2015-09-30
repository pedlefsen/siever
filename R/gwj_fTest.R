#' THIS IS SIMPLY A PLACE HOLDER FOR fTEST WHICH IS NOT YET IMPLEMENTED
#'
#' ITS BODY IS COPY/PASTE OF THE tTEST. TO BE COMPLETED LATER
#'
#' Looks like an internal function to me (Ted).  Not gonna export it unless
#' it turns out to be necessary later.
#' 
#' @param weights A matrix.
#' @param vaccine.indices indices into \code{weight} for vaccine recipients.
#' @param placebo.indices indices into \code{weight} for placebo recipients.
GWJ_fTest <-
  function( weights,vaccine.indices, placebo.indices){
    cat("Not implemented uses the test test as the body for now")
    var.vaccine.vec <- apply(weights[vaccine.indices, ],MARGIN=2, FUN=var, na.rm = T)
    var.placebo.vec <- apply(weights[placebo.indices, ],MARGIN=2, FUN=var, na.rm = T)
    
    obs.num.v.weight <- length(vaccine.indices)
    obs.num.p.weight <- length(placebo.indices)
    obs.num.total <- obs.num.v.weight + obs.num.p.weight
    
    pooled.se.vec <- sqrt((obs.num.v.weight - 1) * var.vaccine.vec/(obs.num.total - 2) + 
                            (obs.num.p.weight - 1) * var.placebo.vec/(obs.num.total - 2))
    
    vaccine.mean.weights <- apply(weights[vaccine.indices, ],MARGIN=2,FUN=mean, na.rm = T)
    placebo.mean.weights <- apply(weights[placebo.indices, ],MARGIN=2,FUN=mean, na.rm = T)
    
    tStats.alt <- (vaccine.mean.weights - placebo.mean.weights)/pooled.se.vec
    
    return(tStats.alt)
  }