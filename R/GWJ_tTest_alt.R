#' ALTERNATE IMPLEMENTATION OF GWJ t-statistic driver
#'
#'NOTE: currently if at a site we have all  placebo and all vaccine constant and equal the t.stat is 0/0 = NaN.
#'In the context of the sieve analysis it might make sense to force that to be p=1 but need double checking
#' Blah.  Blah.  Don't know what to say (Ted).  Looks like a very abstract implementation
#' of the algorithm.  Seems to depend purely on a weight matrix.  Someone else must write this.
#' Not sure if this is an internal function.  It looks like one.  Temporarily \emph{not} gonna
#' export it.
#'
#' @param weight.obj A matrix.
#' @param vaccine.indices indices into \code{weight} for vaccine recipients.
#' @param placebo.indices indices into \code{weight} for placebo recipients.
#' @param perm boolean, when TRUE, df is computed with equal variance assumption
GWJ_tTest.alt <- ### TODO: RENAME NO UNDERSCORES
  function ( weight.obj, vaccine.indices, placebo.indices, perm) {
    equal.var.assumption <- perm
    
    var.vaccine.vec <- apply( weight.obj$weights[ vaccine.indices, , drop = FALSE ], MARGIN=2, FUN=var, na.rm = T );
    var.placebo.vec <- apply( weight.obj$weights[ placebo.indices, , drop = FALSE ], MARGIN=2, FUN=var, na.rm = T );
    
    # obs.num.v.weight <- length(vaccine.indices);  ## TODO: Should be sum( !is.na( weights[ vaccine.indices, , drop = FALSE ] ) )
    # obs.num.p.weight <- length(placebo.indices);  ## TODO: same
    obs.num.v.weight <- apply( ( !is.na( weight.obj$weights[ vaccine.indices, , drop = FALSE ] ) ), MARGIN=2, FUN=sum );
    obs.num.p.weight <- apply( ( !is.na( weight.obj$weights[ placebo.indices, , drop = FALSE ] ) ), MARGIN=2, FUN=sum );
    #obs.num.total <- apply( rbind( obs.num.v.weight, obs.num.p.weight ), MARGIN = 2, FUN = sum );
    obs.num.total <- obs.num.v.weight + obs.num.p.weight
    
    pooled.se.vec <-
      sqrt( ( obs.num.v.weight - 1 ) * var.vaccine.vec/( obs.num.total - 2 ) + 
              ( obs.num.p.weight - 1 ) * var.placebo.vec/( obs.num.total - 2 ) );
    pooled.sample.size <- 1/obs.num.v.weight + 1/obs.num.p.weight;
    
    vaccine.mean.weights <- apply( weight.obj$weights[ vaccine.indices, , drop = FALSE ], MARGIN=2, FUN=mean, na.rm = T );
    placebo.mean.weights <- apply( weight.obj$weights[ placebo.indices, , drop = FALSE ], MARGIN=2, FUN=mean, na.rm = T );
    
    tStats <- ( vaccine.mean.weights - placebo.mean.weights ) / ( pooled.se.vec * sqrt( pooled.sample.size ) )
    
    df <- obs.num.total - 2
    if(!equal.var.assumption){# df based on Welch–Satterthwaite equation for unequal variance
      var.over.n.v <- var.vaccine.vec/obs.num.v.weight
      var.over.n.p <- var.placebo.vec/obs.num.p.weight
      df <- (var.over.n.v + var.over.n.p)^2 /(var.over.n.v^2 /(obs.num.v.weight -1) + var.over.n.p^2 /(obs.num.p.weight - 1) )
    }
    
    return(list( test.stats = tStats, df = df ));
  } # GWJ_tTest.alt (..)
